#include "mosflm/diffraction/SummationIntegrator.hpp"

#include "mosflm/diffraction/LpCorrection.hpp"

#include <cmath>
#include <cstdint>

namespace mosflm::diffraction {
namespace {

inline std::int32_t pixel_at(
    const image::DiffractionImage& image,
    int col,
    int row
) {
    // CBFlib / DiffractionImage layout: row-major with slow = height (ndim1),
    // fast = width (ndim2): index = row * width + col.
    const int w = image.metadata.width;
    const int h = image.metadata.height;
    if (col < 0 || row < 0 || col >= w || row >= h) {
        return 0;
    }
    return image.pixels[static_cast<std::size_t>(row) * static_cast<std::size_t>(w) +
                        static_cast<std::size_t>(col)];
}

}  // namespace

void SummationIntegrator::mm_to_image_pixel(
    const DetectorGeometry& detector,
    double x_mm,
    double y_mm,
    int& col,
    int& row
) {
    double xf = 0.0;
    double yf = 0.0;
    detector.mm_to_pixel(x_mm, y_mm, xf, yf);
    col = static_cast<int>(std::lround(xf + 0.5 * detector.nfast));
    row = static_cast<int>(std::lround(yf + 0.5 * detector.nslow));
}

std::vector<IntegratedReflection> SummationIntegrator::integrate(
    const IntegrationRequest& request
) const {
    const auto& image = request.image;
    const auto& detector = request.detector;
    const auto& predictions = request.reflections;
    const auto& settings = request.settings;

    std::vector<IntegratedReflection> out;
    out.reserve(predictions.size());

    const int half = settings.box_half_px;
    const int frame = settings.bg_frame_px;
    if (half < 0 || frame < 0) {
        return out;
    }

    for (const auto& pred : predictions) {
        IntegratedReflection r;
        r.h = pred.h;
        r.k = pred.k;
        r.l = pred.l;
        r.x_mm = pred.x_mm;
        r.y_mm = pred.y_mm;

        int col0 = 0;
        int row0 = 0;
        mm_to_image_pixel(detector, pred.x_mm, pred.y_mm, col0, row0);

        const int outer = half + frame;
        // Background: pixels in outer square but outside peak square.
        double bg_sum = 0.0;
        int n_bg = 0;
        for (int dr = -outer; dr <= outer; ++dr) {
            for (int dc = -outer; dc <= outer; ++dc) {
                if (std::abs(dr) <= half && std::abs(dc) <= half) {
                    continue;
                }
                const int c = col0 + dc;
                const int rr = row0 + dr;
                if (c < 0 || rr < 0 || c >= image.metadata.width ||
                    rr >= image.metadata.height) {
                    continue;
                }
                bg_sum += static_cast<double>(pixel_at(image, c, rr));
                ++n_bg;
            }
        }

        double peak_sum = 0.0;
        int n_peak = 0;
        for (int dr = -half; dr <= half; ++dr) {
            for (int dc = -half; dc <= half; ++dc) {
                const int c = col0 + dc;
                const int rr = row0 + dr;
                if (c < 0 || rr < 0 || c >= image.metadata.width ||
                    rr >= image.metadata.height) {
                    continue;
                }
                peak_sum += static_cast<double>(pixel_at(image, c, rr));
                ++n_peak;
            }
        }

        r.npix_peak = n_peak;
        r.npix_bg = n_bg;
        if (n_peak == 0) {
            r.ok = false;
            out.push_back(r);
            continue;
        }

        const double bg_mean = (n_bg > 0) ? (bg_sum / static_cast<double>(n_bg)) : 0.0;
        r.intensity = peak_sum - bg_mean * static_cast<double>(n_peak);
        // Simple variance model: peak Poisson + background contribution.
        const double var =
            std::fabs(peak_sum) +
            (n_bg > 0
                 ? (static_cast<double>(n_peak * n_peak) / static_cast<double>(n_bg)) *
                       std::fabs(bg_mean)
                 : 0.0);
        r.sigma = std::sqrt(std::max(var, 1.0));
        if (settings.apply_lp) {
            apply_lp_inplace(detector, r.x_mm, r.y_mm, r.intensity, r.sigma);
        }
        r.ok = true;
        out.push_back(r);
    }
    return out;
}

}  // namespace mosflm::diffraction
