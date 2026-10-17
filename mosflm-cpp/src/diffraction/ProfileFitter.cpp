#include "mosflm/diffraction/ProfileFitter.hpp"

#include "mosflm/diffraction/LpCorrection.hpp"

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

struct FitResult {
    double A{0.0};
    double B{0.0};
    double chi2{0.0};
    double intensity{0.0};
    int npix{0};
    bool ok{false};
};

FitResult fit_fixed_sigma(
    const image::DiffractionImage& image,
    int col0,
    int row0,
    int half,
    double sigma_px
) {
    FitResult out;
    if (!(sigma_px > 1e-6) || half < 0) {
        return out;
    }
    const double inv_2s2 = 1.0 / (2.0 * sigma_px * sigma_px);

    // Linear LS: p ≈ B + A w, w = exp(-r²/(2σ²))
    double sw = 0.0;
    double sw2 = 0.0;
    double sp = 0.0;
    double spw = 0.0;
    int n = 0;
    std::vector<double> ws;
    std::vector<double> ps;
    ws.reserve(static_cast<std::size_t>((2 * half + 1) * (2 * half + 1)));
    ps.reserve(ws.capacity());

    for (int dr = -half; dr <= half; ++dr) {
        for (int dc = -half; dc <= half; ++dc) {
            const int c = col0 + dc;
            const int r = row0 + dr;
            if (c < 0 || r < 0 || c >= image.metadata.width ||
                r >= image.metadata.height) {
                continue;
            }
            const double rr = static_cast<double>(dc * dc + dr * dr);
            const double w = std::exp(-rr * inv_2s2);
            const double p = static_cast<double>(pixel_at(image, c, r));
            ws.push_back(w);
            ps.push_back(p);
            sw += w;
            sw2 += w * w;
            sp += p;
            spw += p * w;
            ++n;
        }
    }
    if (n < 3) {
        return out;
    }

    const double nn = static_cast<double>(n);
    const double det = nn * sw2 - sw * sw;
    if (std::fabs(det) < 1e-12) {
        return out;
    }
    const double A = (nn * spw - sw * sp) / det;
    const double B = (sw2 * sp - sw * spw) / det;

    double chi2 = 0.0;
    double wsum = 0.0;
    for (std::size_t i = 0; i < ws.size(); ++i) {
        const double pred = B + A * ws[i];
        const double d = ps[i] - pred;
        chi2 += d * d;
        wsum += ws[i];
    }
    out.A = A;
    out.B = B;
    out.chi2 = chi2 / nn;
    out.intensity = A * wsum;  // discrete profile integral
    out.npix = n;
    out.ok = true;
    return out;
}

}  // namespace

std::vector<ProfileFittedReflection> ProfileFitter::fit(
    const image::DiffractionImage& image,
    const DetectorGeometry& detector,
    const std::vector<PredictedReflection>& predictions,
    const ProfileFitRequest& request
) const {
    std::vector<ProfileFittedReflection> out;
    out.reserve(predictions.size());
    const int half = request.box_half_px;

    for (const auto& pred : predictions) {
        ProfileFittedReflection r;
        r.h = pred.h;
        r.k = pred.k;
        r.l = pred.l;
        r.x_mm = pred.x_mm;
        r.y_mm = pred.y_mm;

        int col0 = 0;
        int row0 = 0;
        SummationIntegrator::mm_to_image_pixel(
            detector, pred.x_mm, pred.y_mm, col0, row0
        );

        double best_sigma = request.sigma_px;
        FitResult best = fit_fixed_sigma(image, col0, row0, half, best_sigma);

        if (request.refine_sigma && best.ok) {
            // Coarse search around nominal σ.
            const double lo = std::max(0.5, request.sigma_px * 0.5);
            const double hi = request.sigma_px * 2.0;
            for (int i = 0; i <= 10; ++i) {
                const double s = lo + (hi - lo) * (static_cast<double>(i) / 10.0);
                const FitResult trial = fit_fixed_sigma(image, col0, row0, half, s);
                if (trial.ok && trial.chi2 < best.chi2) {
                    best = trial;
                    best_sigma = s;
                }
            }
        }

        if (!best.ok) {
            out.push_back(r);
            continue;
        }
        r.amplitude = best.A;
        r.background = best.B;
        r.sigma_px = best_sigma;
        r.intensity = best.intensity;
        r.chi2 = best.chi2;
        r.npix = best.npix;
        // Rough sigma(I): use sqrt(|I| + n*|B|) style.
        r.sigma = std::sqrt(
            std::max(1.0, std::fabs(best.intensity) + std::fabs(best.B) * best.npix)
        );
        if (request.apply_lp) {
            apply_lp_inplace(detector, r.x_mm, r.y_mm, r.intensity, r.sigma);
        }
        r.ok = true;
        out.push_back(r);
    }
    return out;
}

}  // namespace mosflm::diffraction
