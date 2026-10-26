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

namespace {

/// Gaussian weights as ``scipy.ndimage.gaussian_filter`` (truncate 4σ).
std::vector<double> gaussian_weights(double sigma) {
    const int radius = static_cast<int>(4.0 * sigma + 0.5);
    std::vector<double> w(static_cast<std::size_t>(2 * radius + 1));
    double sum = 0.0;
    for (int k = -radius; k <= radius; ++k) {
        w[static_cast<std::size_t>(k + radius)] = std::exp(-0.5 * k * k / (sigma * sigma));
        sum += w[static_cast<std::size_t>(k + radius)];
    }
    for (auto& v : w) {
        v /= sum;
    }
    return w;
}

/// SciPy's ``reflect`` boundary (d c b a | a b c d | d c b a).
int reflect_index(int i, int n) {
    if (n == 1) {
        return 0;
    }
    while (i < 0 || i >= n) {
        i = i < 0 ? -i - 1 : 2 * n - i - 1;
    }
    return i;
}

/// Separable convolution of a ``w`` × ``h`` image with ``k`` along both axes.
std::vector<double> convolve_2d(const std::vector<double>& in, int w, int h, const std::vector<double>& k) {
    const int radius = static_cast<int>(k.size() / 2);
    std::vector<double> tmp(in.size());
    std::vector<double> out(in.size());
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            double sum = 0.0;
            for (int j = -radius; j <= radius; ++j) {
                sum += k[static_cast<std::size_t>(j + radius)] *
                       in[static_cast<std::size_t>(reflect_index(y + j, h)) * w + x];
            }
            tmp[static_cast<std::size_t>(y) * w + x] = sum;
        }
    }
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            double sum = 0.0;
            for (int j = -radius; j <= radius; ++j) {
                sum += k[static_cast<std::size_t>(j + radius)] *
                       tmp[static_cast<std::size_t>(y) * w + reflect_index(x + j, w)];
            }
            out[static_cast<std::size_t>(y) * w + x] = sum;
        }
    }
    return out;
}

/// NaN-aware Gaussian smoothing (normalised convolution); NaN stays NaN.
std::vector<double> smooth_ignoring_nan(const std::vector<double>& img, int w, int h, double sigma) {
    std::vector<double> value(img.size());
    std::vector<double> weight(img.size());
    for (std::size_t i = 0; i < img.size(); ++i) {
        const bool ok = std::isfinite(img[i]);
        value[i] = ok ? img[i] : 0.0;
        weight[i] = ok ? 1.0 : 0.0;
    }
    const auto k = gaussian_weights(sigma);
    const auto num = convolve_2d(value, w, h, k);
    const auto den = convolve_2d(weight, w, h, k);
    std::vector<double> out(img.size());
    for (std::size_t i = 0; i < img.size(); ++i) {
        out[i] = std::isfinite(img[i]) ? num[i] / den[i] : std::numeric_limits<double>::quiet_NaN();
    }
    return out;
}

/// ``numpy.percentile`` with linear interpolation; reorders ``v``.
double percentile(std::vector<double>& v, double p) {
    std::sort(v.begin(), v.end());
    const double pos = p / 100.0 * static_cast<double>(v.size() - 1);
    const auto lo = static_cast<std::size_t>(std::floor(pos));
    const std::size_t hi = std::min(lo + 1, v.size() - 1);
    return v[lo] + (pos - static_cast<double>(lo)) * (v[hi] - v[lo]);
}

/// Pixels of a 4-connected region of ``on`` within ``steps`` grid steps of
/// a seed (``scipy.ndimage.binary_dilation`` with its default cross,
/// iterated ``steps`` times).
std::vector<std::uint8_t> dilate(const std::vector<std::uint8_t>& seed, int w, int h, int steps) {
    std::vector<int> dist(seed.size(), -1);
    std::vector<std::size_t> queue;
    for (std::size_t i = 0; i < seed.size(); ++i) {
        if (seed[i]) {
            dist[i] = 0;
            queue.push_back(i);
        }
    }
    for (std::size_t q = 0; q < queue.size(); ++q) {
        const std::size_t i = queue[q];
        if (dist[i] >= steps) {
            continue;
        }
        const int x = static_cast<int>(i % static_cast<std::size_t>(w));
        const int y = static_cast<int>(i / static_cast<std::size_t>(w));
        const int nx[4] = {x - 1, x + 1, x, x};
        const int ny[4] = {y, y, y - 1, y + 1};
        for (int n = 0; n < 4; ++n) {
            if (nx[n] < 0 || ny[n] < 0 || nx[n] >= w || ny[n] >= h) {
                continue;
            }
            const std::size_t j = static_cast<std::size_t>(ny[n]) * w + nx[n];
            if (dist[j] < 0) {
                dist[j] = dist[i] + 1;
                queue.push_back(j);
            }
        }
    }
    std::vector<std::uint8_t> out(seed.size());
    for (std::size_t i = 0; i < seed.size(); ++i) {
        out[i] = static_cast<std::uint8_t>(dist[i] >= 0);
    }
    return out;
}

}  // namespace

BeamHalo find_beam_halo(
    const std::vector<float>& counts,
    const std::vector<std::uint8_t>& usable,
    const std::vector<std::uint8_t>& shaded,
    int width,
    int height,
    double beam_col,
    double beam_row,
    double row_scale,
    const DefpixSettings& settings
) {
    BeamHalo halo;
    const std::size_t npix = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    if (counts.size() != npix || usable.size() != npix || shaded.size() != npix) {
        return halo;
    }
    const auto ok = [&](std::size_t i) { return usable[i] && counts[i] >= 0.0f; };
    double level = 0.0;
    std::size_t n_level = 0;
    for (std::size_t i = 0; i < npix; ++i) {
        if (ok(i)) {
            level += counts[i];
            ++n_level;
        }
    }
    if (n_level == 0) {
        return halo;
    }
    level /= static_cast<double>(n_level);

    // Bin until the smoothing kernel (its 2σ disc) holds enough counts.
    const double sigma = settings.halo_smooth_px;
    const double needed =
        std::sqrt(settings.halo_kernel_counts / (4.0 * std::numbers::pi * std::max(level, 1e-9)));
    const int size = needed <= sigma ? 1 : std::min(settings.halo_max_bin_px, static_cast<int>(std::ceil(needed)));
    halo.bin_px = size;
    const double kernel_px = 4.0 * std::numbers::pi * std::pow(std::max(sigma, static_cast<double>(size)), 2);
    const int bw = (width + size - 1) / size;
    const int bh = (height + size - 1) / size;
    const std::size_t nb = static_cast<std::size_t>(bw) * static_cast<std::size_t>(bh);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> img(nb, nan);
    std::vector<std::uint8_t> shadow(nb, 0);
    for (int by = 0; by < bh; ++by) {
        for (int bx = 0; bx < bw; ++bx) {
            double sum = 0.0;
            int n = 0;
            int n_shaded = 0;
            for (int r = by * size; r < std::min((by + 1) * size, height); ++r) {
                for (int c = bx * size; c < std::min((bx + 1) * size, width); ++c) {
                    const std::size_t i = static_cast<std::size_t>(r) * static_cast<std::size_t>(width) +
                                          static_cast<std::size_t>(c);
                    if (ok(i)) {
                        sum += counts[i];
                        ++n;
                    }
                    n_shaded += shaded[i] ? 1 : 0;
                }
            }
            const std::size_t b = static_cast<std::size_t>(by) * bw + bx;
            if (2 * n >= size * size) {
                img[b] = sum / n;
            }
            shadow[b] = static_cast<std::uint8_t>(2 * n_shaded > size * size);
        }
    }
    img = smooth_ignoring_nan(img, bw, bh, std::max(1.0, sigma / size));

    // Bin centres' distance from the beam, in bins.
    std::vector<double> radius(nb);
    int n_rings = 0;
    for (int by = 0; by < bh; ++by) {
        for (int bx = 0; bx < bw; ++bx) {
            const double dc = bx * size + 0.5 * (size - 1) - beam_col;
            const double dr = (by * size + 0.5 * (size - 1) - beam_row) * row_scale;
            const double rb = std::hypot(dc, dr) / size;
            radius[static_cast<std::size_t>(by) * bw + bx] = rb;
            n_rings = std::max(n_rings, static_cast<int>(rb) + 1);
        }
    }
    // Each ring's normal level: the median of its unshaded bins, then the
    // lowest over the window of rings further out (SciPy minimum_filter1d,
    // size window + 1, origin −window/2, reflect).
    std::vector<std::vector<double>> per_ring(static_cast<std::size_t>(n_rings));
    for (std::size_t b = 0; b < nb; ++b) {
        if (std::isfinite(img[b]) && !shadow[b]) {
            per_ring[static_cast<std::size_t>(radius[b])].push_back(img[b]);
        }
    }
    std::vector<double> ring_median(static_cast<std::size_t>(n_rings), std::numeric_limits<double>::infinity());
    for (std::size_t k = 0; k < per_ring.size(); ++k) {
        if (!per_ring[k].empty()) {
            ring_median[k] = percentile(per_ring[k], 50.0);
        }
    }
    const int window = std::max(2, static_cast<int>(settings.halo_window_px) / size);
    const int lo = window / 2 - (window + 1) / 2;
    std::vector<double> normal(static_cast<std::size_t>(n_rings), nan);
    for (int k = 0; k < n_rings; ++k) {
        double m = std::numeric_limits<double>::infinity();
        for (int j = k + lo; j <= k + lo + window; ++j) {
            m = std::min(m, ring_median[static_cast<std::size_t>(reflect_index(j, n_rings))]);
        }
        normal[static_cast<std::size_t>(k)] = std::isfinite(m) ? m : nan;
    }

    // Bright: above factor × normal, and significant in the kernel's counts.
    std::vector<std::uint8_t> bright(nb, 0);
    for (std::size_t b = 0; b < nb; ++b) {
        const double v = img[b];
        const double nl = normal[static_cast<std::size_t>(radius[b])];
        if (!std::isfinite(v) || !std::isfinite(nl) || shadow[b]) {
            continue;
        }
        const double excess = (v - nl) * kernel_px;
        const double noise = std::sqrt(std::max(nl, 0.0) * kernel_px);
        bright[b] = static_cast<std::uint8_t>(v > settings.halo_factor * nl &&
                                              excess > settings.halo_significance * std::max(noise, 1.0));
    }
    // Keep the bright regions that touch the beam or come near the shadow.
    const int reach = std::max(2, static_cast<int>(std::ceil(settings.halo_reach_px / size)));
    const double beam_bins = std::max(1.5, settings.halo_beam_px / size);
    const auto near_shadow = dilate(shadow, bw, bh, reach);
    // 1: a region along the shadow; 2: a region at the beam.
    std::vector<std::uint8_t> keep(nb, 0);
    std::vector<std::size_t> region;
    std::vector<std::uint8_t> seen(nb, 0);
    for (std::size_t s0 = 0; s0 < nb; ++s0) {
        if (!bright[s0] || seen[s0]) {
            continue;
        }
        region.assign(1, s0);
        seen[s0] = 1;
        bool at_beam = false;
        bool by_shadow = false;
        for (std::size_t q = 0; q < region.size(); ++q) {
            const std::size_t i = region[q];
            at_beam = at_beam || radius[i] <= beam_bins;
            by_shadow = by_shadow || near_shadow[i];
            const int x = static_cast<int>(i % static_cast<std::size_t>(bw));
            const int y = static_cast<int>(i / static_cast<std::size_t>(bw));
            const int nx[4] = {x - 1, x + 1, x, x};
            const int ny[4] = {y, y, y - 1, y + 1};
            for (int n = 0; n < 4; ++n) {
                if (nx[n] < 0 || ny[n] < 0 || nx[n] >= bw || ny[n] >= bh) {
                    continue;
                }
                const std::size_t j = static_cast<std::size_t>(ny[n]) * bw + nx[n];
                if (bright[j] && !seen[j]) {
                    seen[j] = 1;
                    region.push_back(j);
                }
            }
        }
        if (at_beam || by_shadow) {
            for (const std::size_t i : region) {
                keep[i] = at_beam ? 2 : 1;
            }
        }
    }

    // Back to pixels. The extent around the beam comes from the regions at
    // the beam; patches along a beamstop arm can lie far out.
    std::vector<double> radii;
    halo.mask.assign(npix, 0);
    for (int r = 0; r < height; ++r) {
        for (int c = 0; c < width; ++c) {
            const std::size_t i = static_cast<std::size_t>(r) * static_cast<std::size_t>(width) +
                                  static_cast<std::size_t>(c);
            const std::uint8_t k = keep[static_cast<std::size_t>(r / size) * bw + c / size];
            if (!usable[i] || k == 0) {
                continue;
            }
            halo.mask[i] = 1;
            ++halo.n_pixels;
            if (k == 2) {
                radii.push_back(std::hypot(c - beam_col, (r - beam_row) * row_scale));
            }
        }
    }
    if (halo.n_pixels == 0) {
        halo.mask.clear();
        return halo;
    }
    halo.found = true;
    if (!radii.empty()) {
        halo.radius_px = *std::max_element(radii.begin(), radii.end());
        halo.median_radius_px = percentile(radii, 50.0);
    }
    return halo;
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
            // Overloaded in every sampled frame: no spot stays saturated across
            // the sweep, so this is a gap or a dead pixel the detector reports
            // at its saturation value (16-bit Eiger gaps read 65535). Masked
            // even without the defect search, or a peak touching it is
            // reported as a profile-fitted overload.
            // Otherwise overloaded in half the frames, or in two and near
            // background in the rest (flickering): a spot has bright frames
            // around its saturated ones.
            if (detect && n_over == nf) {
                set(i, PixelStatus::Defective);
            } else if (settings.find_defects && detect && n_over >= 2 &&
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
        // A halo filling whole rings is its own ring median above, so the
        // shaded test cannot see it. Counts summed over the sampled frames
        // (the background leaves out each pixel's largest value).
        std::vector<float> counts(npix, -1.0f);
        std::vector<std::uint8_t> usable(npix);
        std::vector<std::uint8_t> shaded(npix);
        for (std::size_t i = 0; i < npix; ++i) {
            usable[i] = static_cast<std::uint8_t>(is(i, PixelStatus::Trusted) && sums.n_valid[i] >= 2);
            shaded[i] = static_cast<std::uint8_t>(is(i, PixelStatus::Shaded));
            if (usable[i] && res.background[i] >= 0.0f) {
                counts[i] = res.background[i] * static_cast<float>(sums.n_valid[i] - 1);
            }
        }
        res.halo = find_beam_halo(counts, usable, shaded, w, h, bc, br, py / px, settings);
        if (res.halo.found) {
            for (std::size_t i = 0; i < npix; ++i) {
                if (res.halo.mask[i] && is(i, PixelStatus::Trusted)) {
                    set(i, PixelStatus::Halo);
                }
            }
        }
        if (res.halo.radius_px > 0.0) {
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
