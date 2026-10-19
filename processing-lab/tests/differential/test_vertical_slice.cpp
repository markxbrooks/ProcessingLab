#include "mosflm/crystal/MtzAdapter.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/diffraction/SummationIntegrator.hpp"
#include "mosflm/legacy/LegacyMosflmAdapter.hpp"
#include "test_util.hpp"

#include <cstdint>
#include <filesystem>
#include <vector>

using mosflm::crystal::MtzAdapter;
using mosflm::crystal::SimpleReflection;
using mosflm::diffraction::IntegrationRequest;
using mosflm::diffraction::IntegrationSettings;
using mosflm::diffraction::ReflectionPredictor;
using mosflm::diffraction::SummationIntegrator;
using mosflm::image::DiffractionImage;

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
    constexpr int half = 2;
    constexpr std::int32_t peak = 40;
    for (const auto& p : preds) {
        int col = 0, row = 0;
        SummationIntegrator::mm_to_image_pixel(req.detector, p.x_mm, p.y_mm, col,
                                               row);
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
                             static_cast<std::size_t>(c)] += peak;
            }
        }
    }

    IntegrationSettings settings;
    settings.box_half_px = half;
    settings.bg_frame_px = 1;
    const IntegrationRequest ireq{image, req.detector, preds, settings};
    const auto integrated =
        mosflm::legacy::LegacyMosflmAdapter::integrate_summation(ireq);

    std::vector<SimpleReflection> simple;
    for (const auto& r : integrated) {
        if (!r.ok) {
            continue;
        }
        simple.push_back({r.h, r.k, r.l, static_cast<float>(r.intensity),
                          static_cast<float>(r.sigma)});
    }
    CHECK(!simple.empty());

    const std::filesystem::path mtz =
        std::filesystem::temp_directory_path() / "mosflm_vertical_slice.mtz";
    MtzAdapter::write_simple(mtz, req.crystal, simple);
    CHECK(std::filesystem::exists(mtz));

    const auto back = MtzAdapter::read_simple(mtz);
    CHECK_EQ(back.size(), simple.size());
    if (!back.empty()) {
        CHECK_EQ(back[0].h, simple[0].h);
        CHECK_EQ(back[0].k, simple[0].k);
        CHECK_EQ(back[0].l, simple[0].l);
        CHECK_CLOSE(back[0].intensity, simple[0].intensity, 1e-3);
    }

    std::filesystem::remove(mtz);
    return mosflm_test::finish("vertical_slice");
}
