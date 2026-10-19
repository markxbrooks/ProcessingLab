#include "mosflm/diffraction/DetectorInit.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace mosflm::diffraction {
namespace {

void append_corner(
    const image::DiffractionImage& image,
    int c0,
    int r0,
    int box,
    std::vector<double>& out
) {
    const int w = image.metadata.width;
    const int h = image.metadata.height;
    const int c1 = std::min(w, c0 + box);
    const int r1 = std::min(h, r0 + box);
    for (int r = std::max(0, r0); r < r1; ++r) {
        for (int c = std::max(0, c0); c < c1; ++c) {
            const auto v =
                image.pixels[static_cast<std::size_t>(r) *
                                 static_cast<std::size_t>(w) +
                             static_cast<std::size_t>(c)];
            if (v > 0 && v < 60000) {
                out.push_back(static_cast<double>(v));
            }
        }
    }
}

double median_of(std::vector<double> v) {
    if (v.empty()) {
        return 0.0;
    }
    const std::size_t mid = v.size() / 2U;
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(mid), v.end());
    return v[mid];
}

double mad_of(const std::vector<double>& v, double med) {
    if (v.empty()) {
        return 0.0;
    }
    std::vector<double> absdev;
    absdev.reserve(v.size());
    for (double x : v) {
        absdev.push_back(std::fabs(x - med));
    }
    return median_of(std::move(absdev));
}

}  // namespace

DetectorInitEstimate estimate_blank_gain(
    const image::DiffractionImage& image,
    int corner_box_px
) {
    DetectorInitEstimate out;
    const int w = image.metadata.width;
    const int h = image.metadata.height;
    if (w < 8 || h < 8 || image.pixels.empty()) {
        return out;
    }
    const int box = std::clamp(corner_box_px, 4, std::min(w, h) / 4);
    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(4 * box * box));
    append_corner(image, 0, 0, box, samples);
    append_corner(image, w - box, 0, box, samples);
    append_corner(image, 0, h - box, box, samples);
    append_corner(image, w - box, h - box, box, samples);
    if (samples.size() < 16) {
        return out;
    }

    const double med = median_of(samples);
    const double mad = mad_of(samples, med);
    const double sigma = mad * 1.4826;
    out.blank_adu = med;
    out.bg_sigma_adu = std::max(sigma, 1.0);
    // GAIN ≈ var / mean for Poisson-like background (XDS INIT).
    const double mean = std::max(med - 0.0, 1.0);
    const double var = sigma * sigma;
    out.gain = std::clamp(var / mean, 0.25, 8.0);
    out.ok = true;
    return out;
}

}  // namespace mosflm::diffraction
