#include "mosflm/diffraction/Defpix.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace mosflm::diffraction {
namespace {

std::string upper(std::string s) {
    for (auto& ch : s) {
        ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    }
    return s;
}

/// Per-pixel statistics over the sampled frames.
struct PixelSums {
    std::vector<float> sum;
    std::vector<float> max;
    std::vector<float> min;
    std::vector<std::uint8_t> n_valid;
    std::vector<std::uint8_t> n_overloaded;

    explicit PixelSums(std::size_t n)
        : sum(n, 0.0f),
          max(n, -std::numeric_limits<float>::infinity()),
          min(n, std::numeric_limits<float>::infinity()),
          n_valid(n, 0),
          n_overloaded(n, 0) {}

    void add(const std::vector<std::int32_t>& pixels) {
        for (std::size_t i = 0; i < pixels.size(); ++i) {
            const std::int32_t v = pixels[i];
            if (v < 0) {
                if (v == image::kPixelOverloaded) {
                    ++n_overloaded[i];
                }
                continue;
            }
            const auto f = static_cast<float>(v);
            sum[i] += f;
            max[i] = std::max(max[i], f);
            min[i] = std::min(min[i], f);
            ++n_valid[i];
        }
    }

    /// Mean counts per frame without the pixel's largest value; −1 if none.
    [[nodiscard]] float background(std::size_t i) const {
        const int n = n_valid[i];
        if (n == 0) {
            return -1.0f;
        }
        return n == 1 ? sum[i] : (sum[i] - max[i]) / static_cast<float>(n - 1);
    }
};

double median(std::vector<double> v) {
    if (v.empty()) {
        return 0.0;
    }
    const auto mid = v.begin() + static_cast<std::ptrdiff_t>(v.size() / 2);
    std::nth_element(v.begin(), mid, v.end());
    return *mid;
}

}  // namespace

bool UntrustedRegion::contains(int col, int row) const {
    // XDS pixel numbers are 1-based.
    const double x = static_cast<double>(col) + 1.0;
    const double y = static_cast<double>(row) + 1.0;
    const auto& v = values;
    switch (shape) {
        case Shape::Rectangle:
            return x >= std::min(v[0], v[1]) && x <= std::max(v[0], v[1]) && y >= std::min(v[2], v[3]) &&
                   y <= std::max(v[2], v[3]);
        case Shape::Ellipse: {
            const double ax = 0.5 * std::fabs(v[1] - v[0]);
            const double ay = 0.5 * std::fabs(v[3] - v[2]);
            if (ax <= 0.0 || ay <= 0.0) {
                return false;
            }
            const double dx = (x - 0.5 * (v[0] + v[1])) / ax;
            const double dy = (y - 0.5 * (v[2] + v[3])) / ay;
            return dx * dx + dy * dy <= 1.0;
        }
        case Shape::Quadrilateral: {
            // Crossing number; points on an edge count as inside.
            bool inside = false;
            for (int i = 0, j = 3; i < 4; j = i++) {
                const double xi = v[2 * i];
                const double yi = v[2 * i + 1];
                const double xj = v[2 * j];
                const double yj = v[2 * j + 1];
                const double cross = (xj - xi) * (y - yi) - (yj - yi) * (x - xi);
                if (std::fabs(cross) < 1e-9 && x >= std::min(xi, xj) && x <= std::max(xi, xj) &&
                    y >= std::min(yi, yj) && y <= std::max(yi, yj)) {
                    return true;
                }
                if ((yi > y) != (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) {
                    inside = !inside;
                }
            }
            return inside;
        }
    }
    return false;
}

std::string UntrustedRegion::shape_name() const {
    switch (shape) {
        case Shape::Rectangle:
            return "RECTANGLE";
        case Shape::Ellipse:
            return "ELLIPSE";
        case Shape::Quadrilateral:
            return "QUADRILATERAL";
    }
    return {};
}

UntrustedRegion parse_untrusted_region(const std::vector<std::string>& args) {
    if (args.empty()) {
        throw std::invalid_argument("UNTRUSTED needs RECTANGLE, ELLIPSE or QUADRILATERAL and its values");
    }
    const std::string shape = upper(args[0]);
    UntrustedRegion r;
    std::size_t n = 4;
    if (shape.rfind("RECT", 0) == 0) {
        r.shape = UntrustedRegion::Shape::Rectangle;
    } else if (shape.rfind("ELLI", 0) == 0) {
        r.shape = UntrustedRegion::Shape::Ellipse;
    } else if (shape.rfind("QUAD", 0) == 0) {
        r.shape = UntrustedRegion::Shape::Quadrilateral;
        n = 8;
    } else {
        throw std::invalid_argument("UNTRUSTED: unknown shape '" + args[0] +
                                    "' (RECTANGLE, ELLIPSE or QUADRILATERAL)");
    }
    if (args.size() != n + 1) {
        throw std::invalid_argument("UNTRUSTED " + r.shape_name() + " needs " + std::to_string(n) +
                                    (n == 4 ? " values: x1 x2 y1 y2" : " values: x1 y1 x2 y2 x3 y3 x4 y4"));
    }
    for (std::size_t i = 1; i <= n; ++i) {
        try {
            r.values.push_back(std::stod(args[i]));
        } catch (const std::exception&) {
            throw std::invalid_argument("UNTRUSTED " + r.shape_name() + ": '" + args[i] + "' is not a number");
        }
    }
    return r;
}

bool ResolutionRange::contains(double d_A) const {
    return d_A >= std::min(d1_A, d2_A) && d_A <= std::max(d1_A, d2_A);
}

std::vector<ResolutionRange> ice_ring_ranges() {
    return {{3.93, 3.87},   {3.70, 3.64},   {3.47, 3.41},   {2.70, 2.64},   {2.28, 2.22},
            {2.102, 2.042}, {1.978, 1.918}, {1.948, 1.888}, {1.913, 1.853}, {1.751, 1.691}};
}

std::size_t DefpixResult::count(PixelStatus s) const {
    return static_cast<std::size_t>(std::count(status.begin(), status.end(), static_cast<std::uint8_t>(s)));
}

std::vector<std::uint8_t> DefpixResult::trusted_mask() const {
    std::vector<std::uint8_t> mask(status.size());
    std::transform(status.begin(), status.end(), mask.begin(), [](std::uint8_t s) {
        const auto st = static_cast<PixelStatus>(s);
        return static_cast<std::uint8_t>(st == PixelStatus::Trusted || st == PixelStatus::DetectorMasked);
    });
    return mask;
}

std::size_t DefpixResult::n_rejected() const {
    return status.size() - count(PixelStatus::Trusted) - count(PixelStatus::DetectorMasked);
}

std::vector<std::uint8_t> DefpixResult::spot_search_mask() const {
    std::vector<std::uint8_t> mask(status.size());
    std::transform(status.begin(), status.end(), mask.begin(), [](std::uint8_t s) {
        const auto st = static_cast<PixelStatus>(s);
        return static_cast<std::uint8_t>(st == PixelStatus::Trusted || st == PixelStatus::DetectorMasked ||
                                         st == PixelStatus::ExcludedResolution || st == PixelStatus::Halo);
    });
    return mask;
}

BeamHalo find_beam_halo(
    const std::vector<float>& background,
    const std::vector<std::uint8_t>& usable,
    int width,
    int height,
    double beam_col,
    double beam_row,
    double row_scale,
    int n_frames,
    const DefpixSettings& settings
) {
    BeamHalo halo;
    const std::size_t npix = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    if (background.size() != npix || usable.size() != npix || n_frames < 1) {
        return halo;
    }
    double level = 0.0;
    std::size_t n_level = 0;
    for (std::size_t i = 0; i < npix; ++i) {
        if (usable[i] && background[i] >= 0.0f) {
            level += background[i];
            ++n_level;
        }
    }
    if (n_level == 0) {
        return halo;
    }
    level /= static_cast<double>(n_level);
    // Rings wide enough that one 50 px from the beam holds ~500 counts.
    constexpr double kRingCounts = 500.0;
    constexpr double kRingRadiusPx = 50.0;
    const double frames = std::max(1, n_frames - 1);
    const int ring_px = std::clamp(
        static_cast<int>(std::ceil(kRingCounts / (2.0 * std::numbers::pi * kRingRadiusPx * std::max(level, 1e-9) * frames))), 2, 16
    );
    halo.ring_px = ring_px;

    // The walk may run a window past the search radius, and a few windows on.
    const double reach_px = settings.halo_search_px + 4.0 * settings.halo_window_px;
    const int n_rings = static_cast<int>(reach_px / ring_px) + 1;
    std::vector<std::vector<float>> rings(static_cast<std::size_t>(n_rings));
    const int r0 = std::max(0, static_cast<int>(beam_row - reach_px / std::max(row_scale, 1e-9)) - 1);
    const int r1 = std::min(height - 1, static_cast<int>(beam_row + reach_px / std::max(row_scale, 1e-9)) + 1);
    const int c0 = std::max(0, static_cast<int>(beam_col - reach_px) - 1);
    const int c1 = std::min(width - 1, static_cast<int>(beam_col + reach_px) + 1);
    for (int r = r0; r <= r1; ++r) {
        for (int c = c0; c <= c1; ++c) {
            const std::size_t i = static_cast<std::size_t>(r) * static_cast<std::size_t>(width) +
                                  static_cast<std::size_t>(c);
            if (!usable[i] || background[i] < 0.0f) {
                continue;
            }
            const int k = static_cast<int>(std::hypot(c - beam_col, (r - beam_row) * row_scale) / ring_px);
            if (k < n_rings) {
                rings[static_cast<std::size_t>(k)].push_back(background[i]);
            }
        }
    }
    // Each ring's mean without its brightest 5% (spots left after taking
    // out each pixel's largest frame); too-small rings are left unknown.
    constexpr std::size_t kMinRingPixels = 100;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> profile(static_cast<std::size_t>(n_rings), nan);
    for (std::size_t k = 0; k < rings.size(); ++k) {
        auto& v = rings[k];
        if (v.size() < kMinRingPixels) {
            continue;
        }
        const auto keep = static_cast<std::size_t>(0.95 * static_cast<double>(v.size()));
        std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(keep), v.end());
        double sum = 0.0;
        for (std::size_t j = 0; j < keep; ++j) {
            sum += v[j];
        }
        profile[k] = sum / static_cast<double>(keep);
    }
    // Smooth over three rings.
    std::vector<double> smooth(profile.size(), nan);
    for (std::size_t k = 0; k < profile.size(); ++k) {
        if (std::isnan(profile[k])) {
            continue;
        }
        double sum = 0.0;
        int n = 0;
        for (std::size_t j = k == 0 ? 0 : k - 1; j <= std::min(k + 1, profile.size() - 1); ++j) {
            if (!std::isnan(profile[j])) {
                sum += profile[j];
                ++n;
            }
        }
        smooth[k] = sum / n;
    }
    // Floor: the lowest level within the window further out.
    const int window = std::max(2, static_cast<int>(settings.halo_window_px / ring_px));
    std::vector<double> floor(smooth.size(), nan);
    for (std::size_t k = 0; k < smooth.size(); ++k) {
        double lo = std::numeric_limits<double>::infinity();
        for (std::size_t j = k; j <= std::min(k + static_cast<std::size_t>(window), smooth.size() - 1); ++j) {
            if (!std::isnan(smooth[j])) {
                lo = std::min(lo, smooth[j]);
            }
        }
        floor[k] = std::isfinite(lo) ? lo : nan;
    }
    // Peak within the search radius, then walk out to where the excess over
    // the floor has fallen to the end fraction of the peak's.
    const int search = std::min(n_rings - 1, static_cast<int>(settings.halo_search_px / ring_px));
    int peak = -1;
    for (int k = 0; k <= search; ++k) {
        if (!std::isnan(smooth[static_cast<std::size_t>(k)]) &&
            (peak < 0 || smooth[static_cast<std::size_t>(k)] > smooth[static_cast<std::size_t>(peak)])) {
            peak = k;
        }
    }
    if (peak < 0 || std::isnan(floor[static_cast<std::size_t>(peak)])) {
        return halo;
    }
    const double top = smooth[static_cast<std::size_t>(peak)];
    const double base = floor[static_cast<std::size_t>(peak)];
    halo.peak_radius_px = (peak + 0.5) * ring_px;
    halo.peak_level = top;
    halo.floor_level = base;
    if (!(top >= settings.halo_min_excess * base) || !(top > base)) {
        return halo;
    }
    const double amplitude = top - base;
    for (int k = peak; k < n_rings - window; ++k) {
        const double v = smooth[static_cast<std::size_t>(k)];
        const double f = floor[static_cast<std::size_t>(k)];
        if (!std::isnan(v) && !std::isnan(f) && v - f <= settings.halo_end_fraction * amplitude) {
            halo.found = true;
            halo.radius_px = static_cast<double>(k + 1) * ring_px;
            return halo;
        }
    }
    return halo;  // never flattened within reach: not a halo
}

std::vector<int> defpix_frames(int first, int last, int n) {
    std::vector<int> out;
    const int total = last - first + 1;
    if (total <= 0 || n <= 0) {
        return out;
    }
    if (n >= total) {
        for (int f = first; f <= last; ++f) {
            out.push_back(f);
        }
        return out;
    }
    for (int k = 0; k < n; ++k) {
        const int f = first + static_cast<int>(std::lround(static_cast<double>(k) * (total - 1) / (n - 1)));
        if (out.empty() || out.back() != f) {
            out.push_back(f);
        }
    }
    return out;
}

DefpixResult run_defpix(
    const std::function<image::DiffractionImage(int)>& read_frame,
    const std::vector<int>& frames,
    const DetectorGeometry& detector,
    double beam_x_mm,
    double beam_y_mm,
    const DefpixSettings& settings
) {
    if (frames.empty()) {
        throw std::invalid_argument("DEFPIX: no frames to sample");
    }
    DefpixResult res;
    std::size_t npix = 0;
    PixelSums sums(0);
    for (const int f : frames) {
        const auto image = read_frame(f);
        if (res.frames_used.empty()) {
            res.width = image.metadata.width;
            res.height = image.metadata.height;
            npix = static_cast<std::size_t>(res.width) * static_cast<std::size_t>(res.height);
            sums = PixelSums(npix);
        }
        if (image.metadata.width != res.width || image.metadata.height != res.height ||
            image.pixels.size() != npix) {
            throw std::runtime_error("DEFPIX: frame " + std::to_string(f) + " differs in size from frame " +
                                     std::to_string(frames.front()));
        }
        sums.add(image.pixels);
        res.frames_used.push_back(f);
    }
    const int w = res.width;
    const int h = res.height;
    const int nf = static_cast<int>(res.frames_used.size());
    const bool detect = nf >= 3;
    res.status.assign(npix, static_cast<std::uint8_t>(PixelStatus::Trusted));
    res.background.resize(npix);
    for (std::size_t i = 0; i < npix; ++i) {
        res.background[i] = sums.background(i);
    }

    const double px = std::max(detector.pixel_size_x_mm, 1e-9);
    const double py = std::max(detector.pixel_size_y_mm, 1e-9);
    const double bc = beam_x_mm / px;
    const double br = beam_y_mm / py;
    const double dist = std::max(detector.distance_mm, 1e-6);
    const double wave = std::max(detector.wavelength_A, 1e-6);
    const auto radius_px = [&](int c, int r) { return std::hypot(c - bc, (r - br) * py / px); };
    const auto d_of = [&](int c, int r) {
        const double s = std::sin(0.5 * std::atan2(radius_px(c, r) * px, dist));
        return s < 1e-12 ? 1e6 : wave / (2.0 * s);
    };
    auto set = [&](std::size_t i, PixelStatus s) { res.status[i] = static_cast<std::uint8_t>(s); };
    auto is = [&](std::size_t i, PixelStatus s) { return res.status[i] == static_cast<std::uint8_t>(s); };

    // Detector masks, defects, regions and resolution, in that precedence.
    for (int r = 0; r < h; ++r) {
        for (int c = 0; c < w; ++c) {
            const std::size_t i = static_cast<std::size_t>(r) * static_cast<std::size_t>(w) +
                                  static_cast<std::size_t>(c);
            const int n_over = sums.n_overloaded[i];
            if (sums.n_valid[i] == 0 && n_over == 0) {
                set(i, PixelStatus::DetectorMasked);
                continue;
            }
            if (std::any_of(settings.untrusted.begin(), settings.untrusted.end(),
                            [&](const UntrustedRegion& u) { return u.contains(c, r); })) {
                set(i, PixelStatus::UntrustedRegion);
                continue;
            }
            if (!settings.excluded_resolution.empty()) {
                const double d = d_of(c, r);
                if (std::any_of(settings.excluded_resolution.begin(), settings.excluded_resolution.end(),
                                [&](const ResolutionRange& rr) { return rr.contains(d); })) {
                    set(i, PixelStatus::ExcludedResolution);
                    continue;
                }
            }
            // Overloaded in half the frames, or in two and near background
            // in the rest (flickering): a spot has bright frames around its
            // saturated ones.
            if (settings.find_defects && detect && n_over >= 2 &&
                (2 * n_over >= nf || sums.max[i] < settings.hot_min_counts)) {
                set(i, PixelStatus::Defective);
            } else if (sums.n_valid[i] == 0) {
                set(i, PixelStatus::DetectorMasked);  // overloaded throughout, not judged
            }
        }
    }
    if (!detect || (!settings.find_shadows && !settings.find_defects && !settings.find_halo)) {
        return res;
    }

    // Block size: a typical block should hold shadow_min_counts counts over
    // the sampled frames (fine-pixel, low-background Eiger frames need more
    // pixels than a Pilatus 6M).
    double mean_b = 0.0;
    std::size_t n_b = 0;
    for (std::size_t i = 0; i < npix; ++i) {
        if (is(i, PixelStatus::Trusted) && sums.n_valid[i] >= 2) {
            mean_b += res.background[i];
            ++n_b;
        }
    }
    if (n_b == 0) {
        return res;
    }
    mean_b /= static_cast<double>(n_b);
    const int s = std::clamp(
        static_cast<int>(std::ceil(std::sqrt(settings.shadow_min_counts / std::max(mean_b * nf, 1e-6)))), 8, 32
    );
    res.shadow_block_px = s;

    // Rings of width s about the beam: mean pixel background (for hot pixels).
    const auto ring_of = [&](int c, int r) { return static_cast<int>(radius_px(c, r) / s); };
    const int n_rings = static_cast<int>(std::hypot(w, h) / s) + 2;
    std::vector<double> ring_sum(static_cast<std::size_t>(n_rings), 0.0);
    std::vector<std::size_t> ring_n(static_cast<std::size_t>(n_rings), 0);
    for (int r = 0; r < h; ++r) {
        for (int c = 0; c < w; ++c) {
            const std::size_t i = static_cast<std::size_t>(r) * static_cast<std::size_t>(w) +
                                  static_cast<std::size_t>(c);
            if (is(i, PixelStatus::Trusted) && sums.n_valid[i] >= 2) {
                const auto k = static_cast<std::size_t>(std::min(ring_of(c, r), n_rings - 1));
                ring_sum[k] += res.background[i];
                ++ring_n[k];
            }
        }
    }
    std::vector<double> ring_mean(static_cast<std::size_t>(n_rings), 0.0);
    for (std::size_t k = 0; k < ring_mean.size(); ++k) {
        ring_mean[k] = ring_n[k] > 0 ? ring_sum[k] / static_cast<double>(ring_n[k]) : 0.0;
    }

    if (settings.find_defects) {
        for (int r = 0; r < h; ++r) {
            for (int c = 0; c < w; ++c) {
                const std::size_t i = static_cast<std::size_t>(r) * static_cast<std::size_t>(w) +
                                      static_cast<std::size_t>(c);
                if (!is(i, PixelStatus::Trusted) || sums.n_valid[i] != nf) {
                    continue;
                }
                const double ring = ring_mean[static_cast<std::size_t>(std::min(ring_of(c, r), n_rings - 1))];
                if (sums.min[i] >= std::max(settings.hot_min_counts, settings.hot_factor * ring)) {
                    set(i, PixelStatus::Hot);
                }
            }
        }
    }

    if (settings.find_shadows) {
        // Each block's mean background (pixel values clipped at 10× their
        // ring, so a spot left in one pixel cannot lift the block) is
        // compared with its resolution ring, XDS-style.
        struct Block {
            int c0, r0, ring;
            double value;
            std::size_t n;
        };
        std::vector<Block> blocks;
        for (int r0 = 0; r0 < h; r0 += s) {
            for (int c0 = 0; c0 < w; c0 += s) {
                double sum = 0.0;
                std::size_t n = 0;
                const int cc = std::min(c0 + s / 2, w - 1);
                const int rc = std::min(r0 + s / 2, h - 1);
                const int ring = std::min(ring_of(cc, rc), n_rings - 1);
                const double cap = 10.0 * ring_mean[static_cast<std::size_t>(ring)] + 2.0;
                for (int r = r0; r < std::min(r0 + s, h); ++r) {
                    for (int c = c0; c < std::min(c0 + s, w); ++c) {
                        const std::size_t i = static_cast<std::size_t>(r) * static_cast<std::size_t>(w) +
                                              static_cast<std::size_t>(c);
                        if (is(i, PixelStatus::Trusted) && sums.n_valid[i] >= 2) {
                            sum += std::min(static_cast<double>(res.background[i]), cap);
                            ++n;
                        }
                    }
                }
                if (2 * n >= static_cast<std::size_t>(s) * static_cast<std::size_t>(s)) {
                    blocks.push_back({c0, r0, ring, sum / static_cast<double>(n), n});
                }
            }
        }
        std::vector<std::vector<double>> per_ring(static_cast<std::size_t>(n_rings));
        for (const auto& b : blocks) {
            per_ring[static_cast<std::size_t>(b.ring)].push_back(b.value);
        }
        std::vector<double> ring_median(static_cast<std::size_t>(n_rings), 0.0);
        for (std::size_t k = 0; k < per_ring.size(); ++k) {
            ring_median[k] = median(per_ring[k]);
        }
        // Expected level: the block's own ring, unless the whole ring is
        // shaded (a beamstop filling it): then the highest ring median out
        // to 64 px further, if the ring is below the value range of that.
        const int reach = std::max(1, 64 / s);
        std::vector<double> expected(static_cast<std::size_t>(n_rings), 0.0);
        for (int k = 0; k < n_rings; ++k) {
            double outer = 0.0;
            for (int j = k; j <= std::min(k + reach, n_rings - 1); ++j) {
                outer = std::max(outer, ring_median[static_cast<std::size_t>(j)]);
            }
            const double own = ring_median[static_cast<std::size_t>(k)];
            expected[static_cast<std::size_t>(k)] = 10000.0 * own < settings.value_range_low * outer ? outer : own;
        }
        for (const auto& b : blocks) {
            const double e = expected[static_cast<std::size_t>(b.ring)];
            if (e <= 0.0 || e * nf * static_cast<double>(b.n) < settings.shadow_min_counts) {
                continue;  // too few counts to judge
            }
            const double ratio = 10000.0 * b.value / e;
            if (ratio >= settings.value_range_low && ratio <= settings.value_range_high) {
                continue;
            }
            for (int r = b.r0; r < std::min(b.r0 + s, h); ++r) {
                for (int c = b.c0; c < std::min(b.c0 + s, w); ++c) {
                    const std::size_t i = static_cast<std::size_t>(r) * static_cast<std::size_t>(w) +
                                          static_cast<std::size_t>(c);
                    if (is(i, PixelStatus::Trusted)) {
                        set(i, PixelStatus::Shaded);
                    }
                }
            }
        }
    }

    if (settings.find_halo) {
        // A ring wholly lit by the beam's flare is its own median above, so
        // the shaded test cannot see it: walk out along the radial
        // background instead, from the pixels still trusted.
        std::vector<std::uint8_t> usable(npix);
        for (std::size_t i = 0; i < npix; ++i) {
            usable[i] = static_cast<std::uint8_t>(is(i, PixelStatus::Trusted) && sums.n_valid[i] >= 2);
        }
        res.halo = find_beam_halo(res.background, usable, w, h, bc, br, py / px, nf, settings);
        if (res.halo.found) {
            for (int r = 0; r < h; ++r) {
                for (int c = 0; c < w; ++c) {
                    const std::size_t i = static_cast<std::size_t>(r) * static_cast<std::size_t>(w) +
                                          static_cast<std::size_t>(c);
                    if (is(i, PixelStatus::Trusted) && radius_px(c, r) <= res.halo.radius_px) {
                        set(i, PixelStatus::Halo);
                    }
                }
            }
            const double edge_mm = (res.halo.radius_px + settings.halo_margin_px) * px;
            res.halo_dmax_A = wave / (2.0 * std::sin(0.5 * std::atan2(edge_mm, dist)));
        }
    }
    return res;
}

void write_defpix_cbf(const std::filesystem::path& path, const DefpixResult& result) {
    image::DiffractionImage img;
    img.metadata.width = result.width;
    img.metadata.height = result.height;
    img.pixels.resize(result.status.size());
    const auto nf = static_cast<float>(result.frames_used.size());
    for (std::size_t i = 0; i < img.pixels.size(); ++i) {
        const auto st = static_cast<PixelStatus>(result.status[i]);
        if (st == PixelStatus::DetectorMasked) {
            img.pixels[i] = -1;
        } else if (st != PixelStatus::Trusted) {
            img.pixels[i] = -3;
        } else {
            img.pixels[i] = static_cast<std::int32_t>(std::lround(std::max(0.0f, result.background[i]) * nf));
        }
    }
    // CBF header lines stay under 80 characters (CBFlib folds longer ones).
    std::ostringstream hdr;
    hdr << "# DEFPIX control image written by ProcessingLab\n"
        << "# Trusted: spot-free background summed over " << result.frames_used.size() << " frames\n"
        << "# -3: ignored (untrusted, excluded resolution, shaded, defective, hot, halo)\n"
        << "# -1: masked by the detector\n";
    std::string line = "# Frames";
    for (const int f : result.frames_used) {
        const std::string next = ' ' + std::to_string(f);
        if (line.size() + next.size() > 76) {
            hdr << line << '\n';
            line = "# Frames";
        }
        line += next;
    }
    hdr << line << '\n'
        << "# Untrusted_region " << result.count(PixelStatus::UntrustedRegion) << '\n'
        << "# Excluded_resolution " << result.count(PixelStatus::ExcludedResolution) << '\n'
        << "# Shaded " << result.count(PixelStatus::Shaded) << '\n'
        << "# Defective " << result.count(PixelStatus::Defective) << '\n'
        << "# Hot " << result.count(PixelStatus::Hot) << '\n'
        << "# Halo " << result.count(PixelStatus::Halo) << '\n';
    if (result.halo.found) {
        hdr << "# Halo_radius_px " << std::fixed << std::setprecision(1) << result.halo.radius_px << '\n';
    }
    image::write_minicbf(path, img, hdr.str());
}

}  // namespace mosflm::diffraction
