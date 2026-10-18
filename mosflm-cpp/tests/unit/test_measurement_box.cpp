/**
 * MEAS-style measurement box: peak/bg rim, tilted plane, overlap flag.
 */

#include "mosflm/diffraction/MeasurementBox.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "test_util.hpp"

#include <cstdint>
#include <vector>

using mosflm::DetectorGeometry;
using mosflm::diffraction::IntegrationRequest;
using mosflm::diffraction::IntegrationSettings;
using mosflm::diffraction::MeasIntegrator;
using mosflm::diffraction::MeasurementBoxSettings;
using mosflm::diffraction::PredictedReflection;
using mosflm::diffraction::fit_background_plane;
using mosflm::image::DiffractionImage;

int main() {
    DiffractionImage image;
    image.metadata.width = 64;
    image.metadata.height = 64;
    image.metadata.pixel_size_x_mm = 0.1;
    image.metadata.pixel_size_y_mm = 0.1;
    image.pixels.assign(64u * 64u, 100);

    // Two peaks 4 px apart (will overlap with sep=0.5 mm → 5 px).
    auto plant = [&](int col, int row, std::int32_t v) {
        for (int dr = -1; dr <= 1; ++dr) {
            for (int dc = -1; dc <= 1; ++dc) {
                const int c = col + dc;
                const int r = row + dr;
                image.pixels[static_cast<std::size_t>(r) * 64u +
                             static_cast<std::size_t>(c)] = v;
            }
        }
    };
    plant(20, 20, 500);
    plant(24, 20, 500);

    DetectorGeometry det;
    det.nfast = 64;
    det.nslow = 64;
    det.pixel_size_x_mm = 0.1;
    det.pixel_size_y_mm = 0.1;
    det.distance_mm = 100.0;
    det.wavelength_A = 1.0;

    std::vector<PredictedReflection> preds(2);
    preds[0].h = 1;
    preds[0].k = 0;
    preds[0].l = 0;
    preds[0].x_mm = (20 - 32) * 0.1;
    preds[0].y_mm = (20 - 32) * 0.1;
    preds[1].h = 0;
    preds[1].k = 1;
    preds[1].l = 0;
    preds[1].x_mm = (24 - 32) * 0.1;
    preds[1].y_mm = (20 - 32) * 0.1;

    IntegrationSettings iset;
    iset.box_half_px = 2;
    iset.bg_frame_px = 2;
    MeasurementBoxSettings box;
    box.separation_x_mm = 0.5;  // 5 px
    box.separation_y_mm = 0.5;
    box.peak_half_px = 1;
    box.bg_rim_px = 2;
    box.bg_reject_nsigma = 3.0;
    box.mark_overlaps = true;
    box.integrate_overlaps = true;
    box.use_bg_plane = true;

    const IntegrationRequest ireq{image, det, preds, iset};
    const auto obs = MeasIntegrator{}.integrate(ireq, box);
    CHECK_EQ(static_cast<int>(obs.size()), 2);
    CHECK(obs[0].overlapped);
    CHECK(obs[1].overlapped);
    CHECK(obs[0].ok);
    CHECK(obs[1].ok);
    CHECK(obs[0].intensity > 0.0);
    CHECK(obs[0].sigma > 0.0);
    CHECK(obs[0].npix_bg > 0);

    // Isolated peak should not be overlapped.
    preds.resize(1);
    plant(40, 40, 800);
    preds[0].x_mm = (40 - 32) * 0.1;
    preds[0].y_mm = (40 - 32) * 0.1;
    const IntegrationRequest ireq2{image, det, preds, iset};
    const auto solo = MeasIntegrator{}.integrate(ireq2, box);
    CHECK_EQ(static_cast<int>(solo.size()), 1);
    CHECK(!solo[0].overlapped);
    CHECK(solo[0].ok);
    CHECK(solo[0].intensity > 100.0);

    // Tilted background: OD = 50 + 2*dc + 3*dr outside a clear centre.
    DiffractionImage tilt_img;
    tilt_img.metadata.width = 32;
    tilt_img.metadata.height = 32;
    tilt_img.pixels.assign(32u * 32u, 0);
    for (int r = 0; r < 32; ++r) {
        for (int c = 0; c < 32; ++c) {
            const int dc = c - 16;
            const int dr = r - 16;
            tilt_img.pixels[static_cast<std::size_t>(r) * 32u +
                            static_cast<std::size_t>(c)] =
                static_cast<std::int32_t>(50 + 2 * dc + 3 * dr);
        }
    }
    // Bright peak so rim-only fit still sees the plane.
    for (int dr = -1; dr <= 1; ++dr) {
        for (int dc = -1; dc <= 1; ++dc) {
            tilt_img.pixels[static_cast<std::size_t>(16 + dr) * 32u +
                            static_cast<std::size_t>(16 + dc)] = 500;
        }
    }
    const std::vector<int> empty;
    const auto plane =
        fit_background_plane(tilt_img, 16, 16, 1, 4, empty, empty, 0.0);
    CHECK(plane.ok);
    CHECK_CLOSE(plane.a, 2.0, 1e-6);
    CHECK_CLOSE(plane.b, 3.0, 1e-6);
    CHECK_CLOSE(plane.c, 50.0, 1e-6);

    return mosflm_test::finish("measurement_box");
}
