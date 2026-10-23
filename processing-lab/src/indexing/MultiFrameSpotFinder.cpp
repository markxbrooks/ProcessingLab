#include "mosflm/indexing/MultiFrameSpotFinder.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <utility>

namespace mosflm::indexing {
namespace {

struct Acc {
    double x_sum{0.0};
    double y_sum{0.0};
    double phi_sum{0.0};
    double inten_sum{0.0};
    int n{0};
};

double mid_phi(const Spot& s) {
    return 0.5 * (s.phi_start_deg + s.phi_end_deg);
}

}  // namespace

SpotList find_spots_multiframe(
    int first_frame,
    int last_frame,
    const FrameLoader& load,
    const MultiFrameSpotSettings& settings
) {
    const int n_frames = std::max(1, last_frame - first_frame + 1);
    // Cap per-frame peak search — 50×MAXSPOTS full-detector NMS crawls.
    // Over-sample ×3 so intensity ranking across frames still fills the budget.
    const int per_frame_cap = std::max(
        40, (settings.max_spots * 3 + n_frames - 1) / n_frames
    );
    const int nproc = std::max(1, settings.nproc);

    // Collect per-frame lists so we can stratify by frame (not fragile φ).
    std::vector<std::vector<Spot>> by_frame(
        static_cast<std::size_t>(n_frames)
    );
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
        fs.max_spots = std::min(fs.max_spots, per_frame_cap);
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
    out.header = headers.empty() ? SpotListHeader{} : headers.front();

    // No merge: keep 3D (xy,φ) observations — needed for large-cell DPS.
    if (settings.merge_xy_mm <= 0.0) {
        std::vector<Spot> selected;
        selected.reserve(static_cast<std::size_t>(
            settings.max_spots > 0 ? settings.max_spots : 1500
        ));
        if (settings.max_spots > 0 &&
            static_cast<int>(n_raw) > settings.max_spots) {
            // Equal quota per frame, then fill remainder from global leftovers.
            const int per_frame =
                std::max(1, settings.max_spots / n_frames);
            std::vector<Spot> leftovers;
            leftovers.reserve(n_raw);
            for (auto& frame_spots : by_frame) {
                std::sort(
                    frame_spots.begin(), frame_spots.end(),
                    [](const Spot& a, const Spot& b) {
                        return a.intensity > b.intensity;
                    }
                );
                const std::size_t take = std::min(
                    frame_spots.size(), static_cast<std::size_t>(per_frame)
                );
                for (std::size_t i = 0; i < take; ++i) {
                    selected.push_back(std::move(frame_spots[i]));
                }
                for (std::size_t i = take; i < frame_spots.size(); ++i) {
                    leftovers.push_back(std::move(frame_spots[i]));
                }
            }
            if (static_cast<int>(selected.size()) < settings.max_spots &&
                !leftovers.empty()) {
                std::sort(
                    leftovers.begin(), leftovers.end(),
                    [](const Spot& a, const Spot& b) {
                        return a.intensity > b.intensity;
                    }
                );
                for (auto& s : leftovers) {
                    if (static_cast<int>(selected.size()) >=
                        settings.max_spots) {
                        break;
                    }
                    selected.push_back(std::move(s));
                }
            }
        } else {
            for (auto& frame_spots : by_frame) {
                for (auto& s : frame_spots) {
                    selected.push_back(std::move(s));
                }
            }
            std::sort(
                selected.begin(), selected.end(),
                [](const Spot& a, const Spot& b) {
                    return a.intensity > b.intensity;
                }
            );
            if (settings.max_spots > 0 &&
                static_cast<int>(selected.size()) > settings.max_spots) {
                selected.resize(static_cast<std::size_t>(settings.max_spots));
            }
        }
        out.spots = std::move(selected);
        return out;
    }

    // Merge path: flatten then cluster in (xy,φ).
    std::vector<Spot> raw;
    raw.reserve(n_raw);
    for (auto& frame_spots : by_frame) {
        for (auto& s : frame_spots) {
            raw.push_back(std::move(s));
        }
    }
    std::sort(
        raw.begin(), raw.end(),
        [](const Spot& a, const Spot& b) { return a.intensity > b.intensity; }
    );

    std::vector<Acc> clusters;
    std::vector<char> used(raw.size(), 0);
    const double xy2 = settings.merge_xy_mm * settings.merge_xy_mm;
    for (std::size_t i = 0; i < raw.size(); ++i) {
        if (used[i]) {
            continue;
        }
        Acc acc;
        const auto& seed = raw[i];
        const double seed_phi = mid_phi(seed);
        for (std::size_t j = i; j < raw.size(); ++j) {
            if (used[j]) {
                continue;
            }
            const auto& s = raw[j];
            const double dx = s.x_mm - seed.x_mm;
            const double dy = s.y_mm - seed.y_mm;
            const double dphi = std::fabs(mid_phi(s) - seed_phi);
            if (dx * dx + dy * dy > xy2) {
                continue;
            }
            if (dphi > settings.merge_phi_tol_deg) {
                continue;
            }
            used[j] = 1;
            acc.x_sum += s.x_mm;
            acc.y_sum += s.y_mm;
            acc.phi_sum += mid_phi(s);
            acc.inten_sum += s.intensity;
            ++acc.n;
        }
        if (acc.n >= settings.min_frames) {
            clusters.push_back(acc);
        }
    }

    std::sort(
        clusters.begin(), clusters.end(),
        [](const Acc& a, const Acc& b) { return a.inten_sum > b.inten_sum; }
    );
    if (static_cast<int>(clusters.size()) > settings.max_spots) {
        clusters.resize(static_cast<std::size_t>(settings.max_spots));
    }

    out.spots.reserve(clusters.size());
    for (const auto& c : clusters) {
        Spot s;
        const double inv = 1.0 / static_cast<double>(c.n);
        s.x_mm = c.x_sum * inv;
        s.y_mm = c.y_sum * inv;
        const double phi = c.phi_sum * inv;
        s.phi_start_deg = phi;
        s.phi_end_deg = phi;
        s.intensity = c.inten_sum;
        out.spots.push_back(s);
    }
    return out;
}

}  // namespace mosflm::indexing
