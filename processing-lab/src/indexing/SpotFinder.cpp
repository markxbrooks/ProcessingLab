#include "mosflm/indexing/SpotFinder.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mosflm::indexing {
namespace {

constexpr double kMadScale = 1.4826;
constexpr std::int32_t kMaxValidCounts = 60000;

double percentile_nth(std::vector<double>& v, double frac) {
    if (v.empty()) {
        return 0.0;
    }
    frac = std::clamp(frac, 0.0, 1.0);
    const std::size_t k = static_cast<std::size_t>(
        frac * static_cast<double>(v.size() - 1U)
    );
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(k), v.end());
    return v[k];
}

double median_inplace(std::vector<double>& v) {
    if (v.empty()) {
        return 0.0;
    }
    const std::size_t mid = v.size() / 2U;
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(mid), v.end());
    if ((v.size() % 2U) == 0U && mid > 0U) {
        const double upper = v[mid];
        std::nth_element(
            v.begin(),
            v.begin() + static_cast<std::ptrdiff_t>(mid - 1U),
            v.begin() + static_cast<std::ptrdiff_t>(mid)
        );
        return 0.5 * (v[mid - 1U] + upper);
    }
    return v[mid];
}

bool is_local_maximum(
    const std::vector<float>& work,
    int width,
    int height,
    int x,
    int y
) {
    const float v = work[static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
                         static_cast<std::size_t>(x)];
    for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0) {
                continue;
            }
            const int nx = x + dx;
            const int ny = y + dy;
            if (nx < 0 || ny < 0 || nx >= width || ny >= height) {
                continue;
            }
            const float n =
                work[static_cast<std::size_t>(ny) * static_cast<std::size_t>(width) +
                     static_cast<std::size_t>(nx)];
            if (n > v) {
                return false;
            }
        }
    }
    return true;
}

}  // namespace

SpotList SpotFinder::find(
    const image::DiffractionImage& image,
    const SpotFindSettings& settings
) const {
    const int width = image.metadata.width;
    const int height = image.metadata.height;
    if (width <= 0 || height <= 0) {
        throw std::runtime_error("SpotFinder: invalid image size");
    }
    const std::size_t npix =
        static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    if (image.pixels.size() != npix) {
        throw std::runtime_error("SpotFinder: pixel buffer size mismatch");
    }

    const double pixel_mm = (image.metadata.pixel_size_x_mm > 0.0)
                                ? image.metadata.pixel_size_x_mm
                                : image.metadata.pixel_size_y_mm;
    if (pixel_mm <= 0.0) {
        throw std::runtime_error("SpotFinder: pixel size must be > 0");
    }

    double phi0 = settings.phi_start_deg;
    double phi1 = settings.phi_end_deg;
    if (!(phi1 > phi0) && image.metadata.oscillation_range_deg > 0.0) {
        phi0 = image.metadata.oscillation_start_deg;
        phi1 = phi0 + image.metadata.oscillation_range_deg;
    }
    if (!(phi1 > phi0)) {
        phi1 = phi0;
    }

    std::vector<float> work(npix, 0.0f);
    std::vector<std::uint8_t> valid(npix, 0);
    std::vector<double> sample;
    sample.reserve(npix / 4U);

    const double blank = std::max(0.0, settings.blank_adu);
    const double gain = std::clamp(
        settings.gain > 0.0 ? settings.gain : 1.0, 0.25, 8.0
    );
    for (std::size_t i = 0; i < npix; ++i) {
        const auto v = image.pixels[i];
        if (v > 0 && v < kMaxValidCounts) {
            valid[i] = 1;
            const double adj = static_cast<double>(v) - blank;
            work[i] = static_cast<float>(std::max(0.0, adj));
            sample.push_back(std::max(0.0, adj));
        }
    }
    if (sample.size() < 100U) {
        return {};
    }

    auto sample_bg = sample;
    const double bg_cut = percentile_nth(sample_bg, 0.80);
    std::vector<double> background;
    background.reserve(sample.size());
    for (double v : sample) {
        if (v <= bg_cut) {
            background.push_back(v);
        }
    }
    if (background.size() < 100U) {
        background = sample;
    }

    const double med = median_inplace(background);
    std::vector<double> absdev;
    absdev.reserve(background.size());
    for (double v : background) {
        absdev.push_back(std::fabs(v - med));
    }
    const double mad = median_inplace(absdev) + 1e-6;
    // XDS: var ≈ GAIN · (pixel − blank); inflate MAD σ by √gain.
    const double sigma =
        std::max(kMadScale * mad * std::sqrt(gain), 0.5);
    double cutoff =
        med + static_cast<double>(settings.threshold_sigma) * sigma;
    cutoff = std::max(cutoff, med + 3.0);

    const int border = std::max(0, settings.border);
    const int min_dist = std::max(1, settings.min_distance);
    const int min_dist2 = min_dist * min_dist;
    const int max_spots = std::max(1, settings.max_spots);

    struct Cand {
        int x;
        int y;
        float inten;
    };
    std::vector<Cand> cands;
    cands.reserve(4096);

    for (int y = border; y < height - border; ++y) {
        for (int x = border; x < width - border; ++x) {
            const std::size_t i =
                static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
                static_cast<std::size_t>(x);
            if (!valid[i]) {
                continue;
            }
            if (static_cast<double>(work[i]) < cutoff) {
                continue;
            }
            if (!is_local_maximum(work, width, height, x, y)) {
                continue;
            }
            cands.push_back(Cand{x, y, work[i]});
        }
    }

    std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) {
        return a.inten > b.inten;
    });

    std::vector<Cand> picked;
    picked.reserve(static_cast<std::size_t>(max_spots));
    for (const auto& c : cands) {
        bool ok = true;
        for (const auto& p : picked) {
            const int dx = c.x - p.x;
            const int dy = c.y - p.y;
            if (dx * dx + dy * dy < min_dist2) {
                ok = false;
                break;
            }
        }
        if (!ok) {
            continue;
        }
        picked.push_back(c);
        if (static_cast<int>(picked.size()) >= max_spots) {
            break;
        }
    }

    SpotList out;
    out.header.nfast = width;
    out.header.nslow = height;
    out.header.pixel_mm = pixel_mm;
    out.header.yscale = 1.0;
    out.header.omega_deg = 0.0;
    out.header.invflag = 0;
    out.header.liswung = 1;
    out.header.beam_x_mm = image.metadata.beam_x_mm;
    out.header.beam_y_mm = image.metadata.beam_y_mm;
    out.spots.reserve(picked.size());

    for (const auto& p : picked) {
        Spot s;
        s.x_mm = static_cast<double>(p.x) * pixel_mm;
        s.y_mm = static_cast<double>(p.y) * pixel_mm;
        s.phi_start_deg = phi0;
        s.phi_end_deg = phi1;
        s.intensity = static_cast<double>(p.inten);
        // Use robust background σ (not √I). AUTOINDEX THRESHOLD is min I/σ;
        // √I would cap I/σ at √I and reject typical Eiger peaks at THRESHOLD 20.
        s.sigma = std::max(sigma, 0.5);
        out.spots.push_back(s);
    }
    return out;
}

}  // namespace mosflm::indexing
