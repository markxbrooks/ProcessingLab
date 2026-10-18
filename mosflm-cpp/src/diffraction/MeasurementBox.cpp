#include "mosflm/diffraction/MeasurementBox.hpp"

#include "mosflm/diffraction/LpCorrection.hpp"
#include "mosflm/geometry/LabFrame.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace mosflm::diffraction {
namespace {

inline std::int32_t pixel_at(
    const image::DiffractionImage& image,
    int col,
    int row
) {
    const int w = image.metadata.width;
    const int h = image.metadata.height;
    if (col < 0 || row < 0 || col >= w || row >= h) {
        return 0;
    }
    return image.pixels[static_cast<std::size_t>(row) * static_cast<std::size_t>(w) +
                        static_cast<std::size_t>(col)];
}

void to_pixel(
    const DetectorGeometry& detector,
    const IntegrationSettings& settings,
    double x_mm,
    double y_mm,
    int& col,
    int& row
) {
    if (settings.mosflm_detector_axes) {
        double xdet = 0.0;
        double ydet = 0.0;
        mosflm::cpp_mm_to_mosflm_xdet(
            x_mm,
            y_mm,
            settings.image_beam_x_mm,
            settings.image_beam_y_mm,
            detector.pixel_size_x_mm,
            detector.pixel_size_y_mm,
            xdet,
            ydet
        );
        col = static_cast<int>(std::lround(xdet));
        row = static_cast<int>(std::lround(ydet));
        return;
    }
    SummationIntegrator::mm_to_image_pixel(detector, x_mm, y_mm, col, row);
}

int derive_peak_half(
    const DetectorGeometry& detector,
    const MeasurementBoxSettings& box
) {
    if (box.peak_half_px >= 0) {
        return box.peak_half_px;
    }
    const double px = std::max(detector.pixel_size_x_mm, 1e-9);
    const double py = std::max(detector.pixel_size_y_mm, 1e-9);
    // Peak footprint ≈ half SEPARATION (mm → px), at least 2 px.
    const int hx = static_cast<int>(std::lround(0.5 * box.separation_x_mm / px));
    const int hy = static_cast<int>(std::lround(0.5 * box.separation_y_mm / py));
    return std::max(2, std::min(hx, hy));
}

}  // namespace

std::vector<MeasObservation> MeasIntegrator::integrate(
    const IntegrationRequest& request,
    const MeasurementBoxSettings& box
) const {
    const auto& image = request.image;
    const auto& detector = request.detector;
    const auto& predictions = request.reflections;
    const auto& settings = request.settings;

    std::vector<MeasObservation> out;
    const int n = static_cast<int>(predictions.size());
    if (n == 0) {
        return out;
    }

    const int peak_half = derive_peak_half(detector, box);
    const int rim = std::max(0, box.bg_rim_px);
    const int outer = peak_half + rim;
    if (peak_half < 0) {
        return out;
    }

    std::vector<int> col(static_cast<std::size_t>(n));
    std::vector<int> row(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        to_pixel(
            detector, settings, predictions[static_cast<std::size_t>(i)].x_mm,
            predictions[static_cast<std::size_t>(i)].y_mm, col[static_cast<std::size_t>(i)],
            row[static_cast<std::size_t>(i)]
        );
    }

    // OVERLAP: centres closer than SEPARATION (in pixels).
    const double sep_x_px =
        box.separation_x_mm / std::max(detector.pixel_size_x_mm, 1e-9);
    const double sep_y_px =
        box.separation_y_mm / std::max(detector.pixel_size_y_mm, 1e-9);
    std::vector<char> overlapped(static_cast<std::size_t>(n), 0);
    // Neighbours within ~2*SEPARATION (enough to intrude into the bg rim).
    std::vector<std::vector<int>> neighbours(static_cast<std::size_t>(n));
    if (box.mark_overlaps) {
        const double neigh_x = std::max(sep_x_px, static_cast<double>(2 * outer + 1));
        const double neigh_y = std::max(sep_y_px, static_cast<double>(2 * outer + 1));
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                const double dx = static_cast<double>(
                    col[static_cast<std::size_t>(i)] -
                    col[static_cast<std::size_t>(j)]
                );
                const double dy = static_cast<double>(
                    row[static_cast<std::size_t>(i)] -
                    row[static_cast<std::size_t>(j)]
                );
                if (std::fabs(dx) < sep_x_px && std::fabs(dy) < sep_y_px) {
                    overlapped[static_cast<std::size_t>(i)] = 1;
                    overlapped[static_cast<std::size_t>(j)] = 1;
                }
                if (std::fabs(dx) < neigh_x && std::fabs(dy) < neigh_y) {
                    neighbours[static_cast<std::size_t>(i)].push_back(j);
                    neighbours[static_cast<std::size_t>(j)].push_back(i);
                }
            }
        }
    }

    out.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        MeasObservation r;
        const auto& pred = predictions[static_cast<std::size_t>(i)];
        r.h = pred.h;
        r.k = pred.k;
        r.l = pred.l;
        r.x_mm = pred.x_mm;
        r.y_mm = pred.y_mm;
        r.overlapped = overlapped[static_cast<std::size_t>(i)] != 0;
        r.box_half_px = outer;

        if (r.overlapped && !box.integrate_overlaps) {
            r.ok = false;
            out.push_back(r);
            continue;
        }

        const int c0 = col[static_cast<std::size_t>(i)];
        const int r0 = row[static_cast<std::size_t>(i)];
        const auto& neigh = neighbours[static_cast<std::size_t>(i)];

        // Background candidates: outer square minus peak, minus neighbour peaks.
        std::vector<double> bg_vals;
        bg_vals.reserve(static_cast<std::size_t>((2 * outer + 1) * (2 * outer + 1)));
        for (int dr = -outer; dr <= outer; ++dr) {
            for (int dc = -outer; dc <= outer; ++dc) {
                if (std::abs(dr) <= peak_half && std::abs(dc) <= peak_half) {
                    continue;
                }
                const int c = c0 + dc;
                const int rr = r0 + dr;
                if (c < 0 || rr < 0 || c >= image.metadata.width ||
                    rr >= image.metadata.height) {
                    continue;
                }
                bool in_neighbour = false;
                for (int j : neigh) {
                    if (std::abs(c - col[static_cast<std::size_t>(j)]) <= peak_half &&
                        std::abs(rr - row[static_cast<std::size_t>(j)]) <= peak_half) {
                        in_neighbour = true;
                        break;
                    }
                }
                if (in_neighbour) {
                    continue;
                }
                bg_vals.push_back(static_cast<double>(pixel_at(image, c, rr)));
            }
        }

        // Iterative sigma-clip (BGREJECT stand-in).
        int n_rej = 0;
        if (box.bg_reject_nsigma > 0.0 && bg_vals.size() >= 4) {
            for (int pass = 0; pass < 3; ++pass) {
                double sum = 0.0;
                for (double v : bg_vals) {
                    sum += v;
                }
                const double mean = sum / static_cast<double>(bg_vals.size());
                double var = 0.0;
                for (double v : bg_vals) {
                    const double d = v - mean;
                    var += d * d;
                }
                const double sd =
                    std::sqrt(var / static_cast<double>(bg_vals.size()));
                if (sd < 1e-9) {
                    break;
                }
                const double lim = box.bg_reject_nsigma * sd;
                std::vector<double> kept;
                kept.reserve(bg_vals.size());
                for (double v : bg_vals) {
                    if (std::fabs(v - mean) <= lim) {
                        kept.push_back(v);
                    } else {
                        ++n_rej;
                    }
                }
                if (kept.size() == bg_vals.size() || kept.size() < 4) {
                    bg_vals.swap(kept);
                    break;
                }
                bg_vals.swap(kept);
            }
        }
        r.n_bg_rejected = n_rej;
        r.npix_bg = static_cast<int>(bg_vals.size());

        double peak_sum = 0.0;
        int n_peak = 0;
        for (int dr = -peak_half; dr <= peak_half; ++dr) {
            for (int dc = -peak_half; dc <= peak_half; ++dc) {
                const int c = c0 + dc;
                const int rr = r0 + dr;
                if (c < 0 || rr < 0 || c >= image.metadata.width ||
                    rr >= image.metadata.height) {
                    continue;
                }
                peak_sum += static_cast<double>(pixel_at(image, c, rr));
                ++n_peak;
            }
        }
        r.npix_peak = n_peak;
        if (n_peak == 0) {
            r.ok = false;
            out.push_back(r);
            continue;
        }

        double bg_mean = 0.0;
        if (!bg_vals.empty()) {
            double sum = 0.0;
            for (double v : bg_vals) {
                sum += v;
            }
            bg_mean = sum / static_cast<double>(bg_vals.size());
        }
        r.intensity = peak_sum - bg_mean * static_cast<double>(n_peak);
        const double var =
            std::fabs(peak_sum) +
            (r.npix_bg > 0
                 ? (static_cast<double>(n_peak * n_peak) /
                    static_cast<double>(r.npix_bg)) *
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
