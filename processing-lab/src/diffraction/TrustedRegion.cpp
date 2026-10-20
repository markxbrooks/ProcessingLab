#include "mosflm/diffraction/TrustedRegion.hpp"

#include <algorithm>
#include <cmath>

namespace mosflm::diffraction {

std::vector<std::uint8_t> build_trusted_mask(
    const image::DiffractionImage& image,
    const DetectorGeometry& detector,
    const TrustedRegionSettings& settings
) {
    const int w = image.metadata.width;
    const int h = image.metadata.height;
    std::vector<std::uint8_t> mask(
        static_cast<std::size_t>(std::max(0, w * h)), 1
    );
    if (w <= 0 || h <= 0) {
        return mask;
    }

    const double px = std::max(detector.pixel_size_x_mm, 1e-9);
    const double py = std::max(detector.pixel_size_y_mm, 1e-9);
    double bx = settings.beam_x_mm;
    double by = settings.beam_y_mm;
    if (bx <= 0.0) {
        bx = detector.beam_x_mm > 0.0 ? detector.beam_x_mm
                                     : 0.5 * w * px;
    }
    if (by <= 0.0) {
        by = detector.beam_y_mm > 0.0 ? detector.beam_y_mm
                                     : 0.5 * h * py;
    }
    const double bc = bx / px;
    const double br = by / py;
    const double half_diag =
        0.5 * std::hypot(static_cast<double>(w), static_cast<double>(h));
    const double rmax =
        std::clamp(settings.trusted_radius_frac, 0.1, 1.0) * half_diag;
    const double dist = std::max(detector.distance_mm, 1e-6);
    const double wave = std::max(detector.wavelength_A, 1e-6);

    auto d_of_pixel = [&](int c, int r) -> double {
        const double dx = (static_cast<double>(c) - bc) * px;
        const double dy = (static_cast<double>(r) - br) * py;
        const double rad = std::hypot(dx, dy);
        const double theta = 0.5 * std::atan2(rad, dist);
        const double s2t = std::sin(theta);
        if (s2t < 1e-12) {
            return 1e6;
        }
        return wave / (2.0 * s2t);
    };

    for (int r = 0; r < h; ++r) {
        for (int c = 0; c < w; ++c) {
            const double dc = static_cast<double>(c) - bc;
            const double dr = static_cast<double>(r) - br;
            if (!settings.cover_panel && std::hypot(dc, dr) > rmax) {
                mask[static_cast<std::size_t>(r) * static_cast<std::size_t>(w) +
                     static_cast<std::size_t>(c)] = 0;
                continue;
            }
            if (settings.dmin_A > 0.0 || settings.dmax_A > 0.0) {
                const double d = d_of_pixel(c, r);
                if (settings.dmin_A > 0.0 && d < settings.dmin_A) {
                    mask[static_cast<std::size_t>(r) *
                             static_cast<std::size_t>(w) +
                         static_cast<std::size_t>(c)] = 0;
                    continue;
                }
                if (settings.dmax_A > 0.0 && d > settings.dmax_A) {
                    mask[static_cast<std::size_t>(r) *
                             static_cast<std::size_t>(w) +
                         static_cast<std::size_t>(c)] = 0;
                }
            }
        }
    }
    return mask;
}

}  // namespace mosflm::diffraction
