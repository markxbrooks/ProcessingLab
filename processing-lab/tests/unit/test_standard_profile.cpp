/**
 * Standard-profile (Rossmann) fit: global + multi-region NMASKS tiling.
 */

#include "mosflm/diffraction/StandardProfile.hpp"
#include "test_util.hpp"

#include <cmath>
#include <cstdint>
#include <set>
#include <vector>

using mosflm::DetectorGeometry;
using mosflm::diffraction::IntegrationRequest;
using mosflm::diffraction::IntegrationSettings;
using mosflm::diffraction::PredictedReflection;
using mosflm::diffraction::StandardProfileFitter;
using mosflm::diffraction::StandardProfileSettings;
using mosflm::image::DiffractionImage;

namespace {

void plant_gaussian(
    DiffractionImage& image,
    int c0,
    int r0,
    double sig,
    double amp
) {
    const int w = image.metadata.width;
    const double inv = 1.0 / (2.0 * sig * sig);
    for (int dr = -6; dr <= 6; ++dr) {
        for (int dc = -6; dc <= 6; ++dc) {
            const int c = c0 + dc;
            const int r = r0 + dr;
            if (c < 0 || r < 0 || c >= w || r >= image.metadata.height) {
                continue;
            }
            const double weight =
                std::exp(-static_cast<double>(dc * dc + dr * dr) * inv);
            const std::size_t idx =
                static_cast<std::size_t>(r) * static_cast<std::size_t>(w) +
                static_cast<std::size_t>(c);
            image.pixels[idx] = static_cast<std::int32_t>(
                static_cast<double>(image.pixels[idx]) + amp * weight
            );
        }
    }
}

}  // namespace

int main() {
    DiffractionImage image;
    image.metadata.width = 64;
    image.metadata.height = 64;
    image.metadata.pixel_size_x_mm = 0.1;
    image.metadata.pixel_size_y_mm = 0.1;
    image.pixels.assign(64u * 64u, 50);

    plant_gaussian(image, 32, 32, 1.5, 400.0);

    DetectorGeometry det;
    det.nfast = 64;
    det.nslow = 64;
    det.pixel_size_x_mm = 0.1;
    det.pixel_size_y_mm = 0.1;
    det.distance_mm = 100.0;
    det.wavelength_A = 1.0;

    std::vector<PredictedReflection> preds;
    auto add = [&](int h, int col, int row) {
        PredictedReflection p;
        p.h = h;
        p.k = 0;
        p.l = 0;
        p.x_mm = (col - 32) * 0.1;
        p.y_mm = (row - 32) * 0.1;
        preds.push_back(p);
    };
    for (int i = 0; i < 10; ++i) {
        add(i + 1, 32, 32);
    }
    add(100, 10, 10);
    add(101, 50, 50);

    IntegrationSettings iset;
    iset.apply_lp = false;
    StandardProfileSettings ps;
    ps.box.separation_x_mm = 1.0;
    ps.box.separation_y_mm = 1.0;
    ps.box.peak_half_px = 3;
    ps.box.bg_rim_px = 2;
    ps.box.mark_overlaps = false;
    ps.learn_min_isigi = 3.0;
    ps.min_learners = 5;
    ps.prefer_profile = true;
    ps.n_bins_x = 1;
    ps.n_bins_y = 1;

    const IntegrationRequest ireq{image, det, preds, iset};
    const auto obs = StandardProfileFitter{}.fit(ireq, ps);
    CHECK(obs.size() >= 10);
    CHECK(obs[0].ok);
    CHECK(obs[0].profile_ok);
    CHECK(obs[0].profile_intensity > 100.0);
    CHECK(obs[0].sigma > 0.0);
    CHECK(obs[0].pkratio < 2.0);
    CHECK(obs[0].profile_gof > 0.0);
    CHECK(obs[0].profile_gof < 5.0);
    CHECK_EQ(obs[0].profile_region, 0);

    // Multi-region: two differently shaped peaks in opposite corners.
    DiffractionImage image2;
    image2.metadata = image.metadata;
    image2.pixels.assign(64u * 64u, 50);
    plant_gaussian(image2, 12, 12, 1.2, 500.0);  // sharp
    plant_gaussian(image2, 52, 52, 2.8, 500.0);  // broad

    std::vector<PredictedReflection> preds2;
    auto add2 = [&](int h, int col, int row) {
        PredictedReflection p;
        p.h = h;
        p.k = 0;
        p.l = 0;
        p.x_mm = (col - 32) * 0.1;
        p.y_mm = (row - 32) * 0.1;
        preds2.push_back(p);
    };
    for (int i = 0; i < 10; ++i) {
        add2(i + 1, 12, 12);
    }
    for (int i = 0; i < 10; ++i) {
        add2(i + 101, 52, 52);
    }

    StandardProfileSettings ps2 = ps;
    ps2.n_bins_x = 2;
    ps2.n_bins_y = 2;
    ps2.min_learners = 5;
    const IntegrationRequest ireq2{image2, det, preds2, iset};
    const auto obs2 = StandardProfileFitter{}.fit(ireq2, ps2);
    CHECK_EQ(obs2.size(), 20u);

    std::set<int> regions;
    int n_ok = 0;
    for (const auto& o : obs2) {
        regions.insert(o.profile_region);
        if (o.profile_ok) {
            ++n_ok;
        }
    }
    CHECK(n_ok >= 16);
    // Sharp (first 10) and broad (next 10) land in different 2×2 tiles.
    CHECK(regions.size() >= 2u);
    CHECK(obs2[0].profile_region != obs2[10].profile_region);

    return mosflm_test::finish("standard_profile");
}
