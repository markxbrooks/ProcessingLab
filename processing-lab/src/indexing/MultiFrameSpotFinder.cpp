#include "mosflm/indexing/MultiFrameSpotFinder.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mosflm::indexing {
namespace {

// Per-frame search budget. These are heuristics, not physics: a reflection
// typically shows on a handful of frames, and the finder over-samples so the
// brightest-first cut after merging still fills ``max_spots``. The floor keeps
// short ranges useful; on long ranges it can exceed the total budget (raw
// spots ≤ floor × frames), which costs time but not correctness, since the
// output is capped after merging.
constexpr std::int64_t kFramesPerReflection = 5;
constexpr std::int64_t kOversample = 3;
constexpr std::int64_t kMinPerFrameSpots = 40;
/// A track's next spot must lie within this many radii of its running
/// centroid, so a chain of small steps cannot run on indefinitely. Generous:
/// a reflection near the rotation axis drifts across the frames it spans, and
/// integer peak pixels jitter by ±1 px. One radius split 2.5× more real
/// reflections on 76c1b (115 near-duplicate spots against 46) and raised the
/// indexing φ residual from 0.45° to 0.56°; four radii matched no anchor.
constexpr double kCentroidAnchor = 4.0;
/// Per-frame quota when ``max_spots`` is unlimited and not merging.
constexpr std::size_t kDefaultReserve = 1500;

double mid_phi(const Spot& s) {
    return 0.5 * (s.phi_start_deg + s.phi_end_deg);
}

/// Brightest first; ties by position and φ, so the order never depends on
/// how the frames were scheduled.
bool brighter(const Spot& a, const Spot& b) {
    if (a.intensity != b.intensity) {
        return a.intensity > b.intensity;
    }
    if (a.x_mm != b.x_mm) {
        return a.x_mm < b.x_mm;
    }
    if (a.y_mm != b.y_mm) {
        return a.y_mm < b.y_mm;
    }
    return mid_phi(a) < mid_phi(b);
}

void validate(int first_frame, int last_frame, const MultiFrameSpotSettings& settings) {
    if (last_frame < first_frame) {
        throw std::invalid_argument(
            "find_spots_multiframe: last frame " + std::to_string(last_frame) + " is before first frame " +
            std::to_string(first_frame)
        );
    }
    if (!std::isfinite(settings.merge_xy_mm) || settings.merge_xy_mm < 0.0) {
        throw std::invalid_argument("find_spots_multiframe: merge_xy_mm must be finite and >= 0");
    }
    if (settings.min_frames < 1) {
        throw std::invalid_argument("find_spots_multiframe: min_frames must be >= 1");
    }
}

/// The per-frame cap for ``n_frames`` frames (64-bit, so no overflow).
int per_frame_cap(const MultiFrameSpotSettings& settings, std::int64_t n_frames, bool merging) {
    if (settings.max_spots <= 0) {
        return settings.per_frame.max_spots;  // no total budget
    }
    const std::int64_t factor = merging ? kFramesPerReflection : 1;
    const std::int64_t want =
        (static_cast<std::int64_t>(settings.max_spots) * kOversample * factor + n_frames - 1) / n_frames;
    return static_cast<int>(std::clamp<std::int64_t>(want, kMinPerFrameSpots, INT_MAX));
}

/// A reflection's track: the spots it has taken, one per frame.
struct Track {
    double w_sum{0.0};
    double x_sum{0.0};
    double y_sum{0.0};
    double phi_sum{0.0};
    Spot brightest;
    Spot last;
    std::size_t last_frame{0};
    int n_frames{0};

    void add(const Spot& sp, std::size_t frame) {
        const double w = std::max(sp.intensity, 1e-6);
        w_sum += w;
        x_sum += w * sp.x_mm;
        y_sum += w * sp.y_mm;
        phi_sum += w * mid_phi(sp);
        if (n_frames == 0 || brighter(sp, brightest)) {
            brightest = sp;
        }
        last = sp;
        last_frame = frame;
        ++n_frames;
    }
    [[nodiscard]] double cx() const { return x_sum / w_sum; }
    [[nodiscard]] double cy() const { return y_sum / w_sum; }

    /// The merged spot: weighted position and φ, the brightest frame's
    /// intensity and σ.
    [[nodiscard]] Spot merged() const {
        Spot s;
        s.x_mm = cx();
        s.y_mm = cy();
        const double phi = phi_sum / w_sum;
        s.phi_start_deg = phi;
        s.phi_end_deg = phi;
        s.intensity = brightest.intensity;
        s.sigma = brightest.sigma;
        return s;
    }
};

/// Build tracks frame by frame (see :func:`find_spots_multiframe`).
std::vector<Track> associate(const std::vector<std::vector<Spot>>& by_frame, double radius_mm) {
    const double r2 = radius_mm * radius_mm;
    const auto near = [](double x0, double y0, double x1, double y1, double limit2) {
        const double dx = x1 - x0;
        const double dy = y1 - y0;
        return dx * dx + dy * dy <= limit2;
    };
    const double anchor2 = kCentroidAnchor * kCentroidAnchor * r2;
    std::vector<Track> tracks;
    std::vector<std::size_t> open;  // tracks that took a spot on the previous frame
    for (std::size_t f = 0; f < by_frame.size(); ++f) {
        const auto& spots = by_frame[f];
        // Candidate pairs: within the radius of the track's last spot and the
        // anchor of its centroid. Count them both ways; only sole pairs join.
        std::vector<int> per_track(open.size(), 0);
        std::vector<int> per_spot(spots.size(), 0);
        std::vector<std::size_t> match_of_spot(spots.size(), SIZE_MAX);
        for (std::size_t t = 0; t < open.size(); ++t) {
            const Track& tr = tracks[open[t]];
            for (std::size_t s = 0; s < spots.size(); ++s) {
                if (near(tr.last.x_mm, tr.last.y_mm, spots[s].x_mm, spots[s].y_mm, r2) &&
                    near(tr.cx(), tr.cy(), spots[s].x_mm, spots[s].y_mm, anchor2)) {
                    ++per_track[t];
                    ++per_spot[s];
                    match_of_spot[s] = t;
                }
            }
        }
        std::vector<std::size_t> next_open;
        for (std::size_t s = 0; s < spots.size(); ++s) {
            const std::size_t t = match_of_spot[s];
            if (per_spot[s] == 1 && per_track[t] == 1) {
                tracks[open[t]].add(spots[s], f);
                next_open.push_back(open[t]);
            } else {
                // Unmatched, or ambiguous: a new track rather than a guess.
                Track tr;
                tr.add(spots[s], f);
                tracks.push_back(tr);
                next_open.push_back(tracks.size() - 1);
            }
        }
        open = std::move(next_open);
    }
    return tracks;
}

/// Brightest ``max_spots`` (all if ``≤ 0``) of ``spots``, in order.
void keep_brightest(std::vector<Spot>& spots, int max_spots) {
    std::sort(spots.begin(), spots.end(), brighter);
    if (max_spots > 0 && spots.size() > static_cast<std::size_t>(max_spots)) {
        spots.resize(static_cast<std::size_t>(max_spots));
    }
}

/// Per-frame spots without merging: an equal quota per frame, then the
/// brightest of the rest, so the brightest frames do not starve later φ.
std::vector<Spot> stratified(std::vector<std::vector<Spot>>& by_frame, std::size_t n_raw, int max_spots) {
    std::vector<Spot> selected;
    if (max_spots <= 0 || n_raw <= static_cast<std::size_t>(max_spots)) {
        selected.reserve(n_raw);
        for (auto& frame_spots : by_frame) {
            std::move(frame_spots.begin(), frame_spots.end(), std::back_inserter(selected));
        }
        keep_brightest(selected, max_spots);
        return selected;
    }
    selected.reserve(std::max(static_cast<std::size_t>(max_spots), kDefaultReserve));
    const std::size_t quota = std::max<std::size_t>(1, static_cast<std::size_t>(max_spots) / by_frame.size());
    std::vector<Spot> leftovers;
    leftovers.reserve(n_raw);
    for (auto& frame_spots : by_frame) {
        std::sort(frame_spots.begin(), frame_spots.end(), brighter);
        const std::size_t take = std::min(frame_spots.size(), quota);
        std::move(frame_spots.begin(), frame_spots.begin() + static_cast<std::ptrdiff_t>(take),
                  std::back_inserter(selected));
        std::move(frame_spots.begin() + static_cast<std::ptrdiff_t>(take), frame_spots.end(),
                  std::back_inserter(leftovers));
    }
    std::sort(leftovers.begin(), leftovers.end(), brighter);
    for (auto& s : leftovers) {
        if (selected.size() >= static_cast<std::size_t>(max_spots)) {
            break;
        }
        selected.push_back(std::move(s));
    }
    return selected;
}

}  // namespace

SpotList find_spots_multiframe(
    int first_frame,
    int last_frame,
    const FrameLoader& load,
    const MultiFrameSpotSettings& settings
) {
    validate(first_frame, last_frame, settings);
    const std::int64_t n_frames = static_cast<std::int64_t>(last_frame) - first_frame + 1;
    const bool merging = settings.merge_xy_mm > 0.0;
    // Cap per-frame peak search — 50×MAXSPOTS full-detector NMS crawls.
    const int cap = per_frame_cap(settings, n_frames, merging);
    const int nproc = std::max(1, settings.nproc);

    // Collect per-frame lists so we can stratify by frame (not fragile φ).
    std::vector<std::vector<Spot>> by_frame(static_cast<std::size_t>(n_frames));
    std::vector<SpotListHeader> headers(static_cast<std::size_t>(n_frames));
    // An exception must not leave an OpenMP region (std::terminate): keep
    // the first and rethrow it after the loop.
    std::exception_ptr error;

#if defined(PROCESSINGLAB_HAS_OPENMP) && PROCESSINGLAB_HAS_OPENMP
#pragma omp parallel for num_threads(nproc) schedule(dynamic, 1) \
    if (nproc > 1 && n_frames > 1)
#endif
    for (int f = first_frame; f <= last_frame; ++f) {
      try {
        const auto image = load(f);
        SpotFindSettings fs = settings.per_frame;
        fs.max_spots = std::min(fs.max_spots, cap);
        if (fs.phi_end_deg <= fs.phi_start_deg &&
            image.metadata.oscillation_range_deg > 0.0) {
            fs.phi_start_deg = image.metadata.oscillation_start_deg;
            fs.phi_end_deg =
                fs.phi_start_deg + image.metadata.oscillation_range_deg;
        }
        // Thread-local finder (find is const, but avoid any shared mutable).
        auto list = SpotFinder{}.find(image, fs);
        const std::size_t idx =
            static_cast<std::size_t>(f - first_frame);
        headers[idx] = list.header;
        by_frame[idx] = std::move(list.spots);
      } catch (...) {
#if defined(PROCESSINGLAB_HAS_OPENMP) && PROCESSINGLAB_HAS_OPENMP
#pragma omp critical(multiframe_spot_error)
#endif
        if (!error) {
            error = std::current_exception();
        }
      }
    }
    if (error) {
        std::rethrow_exception(error);
    }

    std::size_t n_raw = 0;
    for (const auto& spots : by_frame) {
        n_raw += spots.size();
    }

    SpotList out;
    // Detector, pixel size and beam: the same for every frame of the series.
    out.header = headers.front();

    // No merge: keep 3D (xy,φ) observations.
    if (!merging) {
        out.spots = stratified(by_frame, n_raw, settings.max_spots);
        return out;
    }

    const auto tracks = associate(by_frame, settings.merge_xy_mm);
    out.spots.reserve(tracks.size());
    for (const auto& tr : tracks) {
        if (tr.n_frames >= settings.min_frames) {
            out.spots.push_back(tr.merged());
        }
    }
    keep_brightest(out.spots, settings.max_spots);
    return out;
}

}  // namespace mosflm::indexing
