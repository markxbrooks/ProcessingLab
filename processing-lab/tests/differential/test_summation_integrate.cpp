#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/diffraction/SummationIntegrator.hpp"
#include "mosflm/image/DiffractionImage.hpp"
#include "mosflm/legacy/LegacyMosflmAdapter.hpp"
#include "test_util.hpp"

#include <cmath>
#include <cstdint>
#include <vector>

using mosflm::diffraction::IntegrationRequest;
using mosflm::diffraction::IntegrationSettings;
using mosflm::diffraction::PredictedReflection;
using mosflm::diffraction::ReflectionPredictor;
using mosflm::diffraction::SummationIntegrator;
using mosflm::image::DiffractionImage;

namespace {

void plant_box(DiffractionImage& image, int col, int row, int half,
               std::int32_t value) {
    for (int dr = -half; dr <= half; ++dr) {
        for (int dc = -half; dc <= half; ++dc) {
            const int c = col + dc;
            const int r = row + dr;
            if (c < 0 || r < 0 || c >= image.metadata.width ||
                r >= image.metadata.height) {
                continue;
            }
            image.pixels[static_cast<std::size_t>(r) *
                             static_cast<std::size_t>(image.metadata.width) +
                         static_cast<std::size_t>(c)] += value;
        }
    }
}

}  // namespace

int main() {
    const auto req = ReflectionPredictor::analytical_still_fixture();
    const auto preds =
        mosflm::legacy::LegacyMosflmAdapter::predict_still(req).reflections;

    DiffractionImage image;
    image.metadata.width = req.detector.nfast;
    image.metadata.height = req.detector.nslow;
    image.pixels.assign(
        static_cast<std::size_t>(image.metadata.width) *
            static_cast<std::size_t>(image.metadata.height),
        0
    );

    constexpr int half = 2;
    constexpr std::int32_t peak_100 = 50;
    constexpr std::int32_t peak_010 = 80;

    std::vector<PredictedReflection> two;
    for (const auto& p : preds) {
        if (p.h == 1 && p.k == 0 && p.l == 0) {
            two.push_back(p);
            int col = 0, row = 0;
            SummationIntegrator::mm_to_image_pixel(req.detector, p.x_mm, p.y_mm,
                                                   col, row);
            plant_box(image, col, row, half, peak_100);
        } else if (p.h == 0 && p.k == 1 && p.l == 0) {
            two.push_back(p);
            int col = 0, row = 0;
            SummationIntegrator::mm_to_image_pixel(req.detector, p.x_mm, p.y_mm,
                                                   col, row);
            plant_box(image, col, row, half, peak_010);
        }
    }
    CHECK_EQ(two.size(), 2u);

    IntegrationSettings settings;
    settings.box_half_px = half;
    settings.bg_frame_px = 1;
    const IntegrationRequest ireq{image, req.detector, two, settings};
    const auto integrated =
        mosflm::legacy::LegacyMosflmAdapter::integrate_summation(ireq);
    CHECK_EQ(integrated.size(), 2u);

    int found = 0;
    for (const auto& r : integrated) {
        CHECK(r.ok);
        CHECK_EQ(r.npix_peak, 25);
        if (r.h == 1 && r.k == 0 && r.l == 0) {
            CHECK_CLOSE(r.intensity, 25.0 * peak_100, 1e-9);
            CHECK_CLOSE(r.sigma, std::sqrt(25.0 * peak_100), 1e-9);
            ++found;
        } else if (r.h == 0 && r.k == 1 && r.l == 0) {
            CHECK_CLOSE(r.intensity, 25.0 * peak_010, 1e-9);
            CHECK_CLOSE(r.sigma, std::sqrt(25.0 * peak_010), 1e-9);
            ++found;
        }
    }
    CHECK_EQ(found, 2);

    return mosflm_test::finish("summation_integrate");
}
