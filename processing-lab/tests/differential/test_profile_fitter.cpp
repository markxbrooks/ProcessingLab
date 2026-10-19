#include "mosflm/diffraction/ProfileFitter.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/diffraction/SummationIntegrator.hpp"
#include "mosflm/legacy/LegacyMosflmAdapter.hpp"
#include "test_util.hpp"

#include <cmath>
#include <cstdint>
#include <vector>

using mosflm::diffraction::ProfileFitRequest;
using mosflm::diffraction::ProfileFitSettings;
using mosflm::diffraction::ProfileFitter;
using mosflm::diffraction::ReflectionPredictor;
using mosflm::diffraction::SummationIntegrator;
using mosflm::image::DiffractionImage;

namespace {

void plant_gaussian(DiffractionImage& image, int col0, int row0, double amp,
                    double bg, double sigma_px, int half) {
    const double inv_2s2 = 1.0 / (2.0 * sigma_px * sigma_px);
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
            image.pixels[static_cast<std::size_t>(r) *
                             static_cast<std::size_t>(image.metadata.width) +
                         static_cast<std::size_t>(c)] =
                static_cast<std::int32_t>(std::lround(bg + amp * w));
        }
    }
}

}  // namespace

int main() {
    auto req = ReflectionPredictor::analytical_still_fixture();
    req.settings.dmin_A = 8.0;
    const auto preds =
        mosflm::legacy::LegacyMosflmAdapter::predict_still(req).reflections;
    CHECK(!preds.empty());

    DiffractionImage image;
    image.metadata.width = req.detector.nfast;
    image.metadata.height = req.detector.nslow;
    image.pixels.assign(
        static_cast<std::size_t>(image.metadata.width) *
            static_cast<std::size_t>(image.metadata.height),
        0
    );

    constexpr int half = 5;
    constexpr double amp = 80.0;
    constexpr double bg = 5.0;
    constexpr double sig = 1.5;
    for (const auto& p : preds) {
        int col = 0, row = 0;
        SummationIntegrator::mm_to_image_pixel(req.detector, p.x_mm, p.y_mm, col,
                                               row);
        plant_gaussian(image, col, row, amp, bg, sig, half);
    }

    ProfileFitSettings settings;
    settings.box_half_px = half;
    settings.sigma_px = sig;
    settings.refine_sigma = true;
    settings.apply_lp = false;
    const ProfileFitRequest freq{image, req.detector, preds, settings};
    const auto fitted = ProfileFitter{}.fit(freq);

    CHECK_EQ(fitted.size(), preds.size());

    int n_ok = 0;
    for (const auto& r : fitted) {
        if (!r.ok) {
            continue;
        }
        ++n_ok;
        CHECK_EQ(r.npix, (2 * half + 1) * (2 * half + 1));
        CHECK_CLOSE(r.amplitude, amp, 3.0);
        CHECK_CLOSE(r.background, bg, 3.0);
        CHECK_CLOSE(r.sigma_px, sig, 0.5);
        CHECK_MSG(r.chi2 < 1.0, "chi2 = " + std::to_string(r.chi2));
        CHECK(r.intensity > 0.0);
    }
    CHECK(n_ok > 0);

    // apply_lp scales the fitted intensity.
    ProfileFitSettings lp_settings = settings;
    lp_settings.refine_sigma = false;
    lp_settings.apply_lp = true;
    const ProfileFitRequest lreq{image, req.detector, preds, lp_settings};
    const auto lp_fitted = ProfileFitter{}.fit(lreq);
    CHECK_EQ(lp_fitted.size(), fitted.size());

    ProfileFitSettings base_settings = settings;
    base_settings.refine_sigma = false;
    base_settings.apply_lp = false;
    const ProfileFitRequest breq{image, req.detector, preds, base_settings};
    const auto base_fitted = ProfileFitter{}.fit(breq);

    CHECK_EQ(lp_fitted.size(), base_fitted.size());
    int scaled = 0;
    for (std::size_t i = 0; i < lp_fitted.size(); ++i) {
        if (!lp_fitted[i].ok || !base_fitted[i].ok) {
            continue;
        }
        if (std::fabs(lp_fitted[i].intensity - base_fitted[i].intensity) > 1e-6) {
            ++scaled;
        }
    }
    CHECK(scaled > 0);

    return mosflm_test::finish("profile_fitter");
}
