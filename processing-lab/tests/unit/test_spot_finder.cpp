/**
 * DiffractionLab-style SpotFinder + SpotList::write_spt round-trip.
 */

#include "mosflm/indexing/SpotFinder.hpp"
#include "mosflm/indexing/SpotList.hpp"
#include "test_util.hpp"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <vector>

using mosflm::image::DiffractionImage;
using mosflm::indexing::SpotFindSettings;
using mosflm::indexing::SpotFinder;
using mosflm::indexing::SpotList;

namespace {

DiffractionImage make_synthetic_frame() {
    DiffractionImage img;
    img.metadata.width = 64;
    img.metadata.height = 64;
    img.metadata.pixel_size_x_mm = 0.1;
    img.metadata.pixel_size_y_mm = 0.1;
    img.metadata.beam_x_mm = 3.2;
    img.metadata.beam_y_mm = 3.2;
    img.metadata.oscillation_start_deg = 0.0;
    img.metadata.oscillation_range_deg = 0.5;
    img.pixels.assign(64U * 64U, 10);
    // Three strong peaks well separated.
    const int peaks[][2] = {{16, 16}, {40, 20}, {30, 48}};
    for (const auto& p : peaks) {
        img.pixels[static_cast<std::size_t>(p[1]) * 64U +
                   static_cast<std::size_t>(p[0])] = 5000;
    }
    return img;
}

}  // namespace

int main() {
    const auto img = make_synthetic_frame();
    SpotFindSettings settings;
    settings.threshold_sigma = 5.0f;
    settings.min_distance = 4;
    settings.max_spots = 20;
    settings.border = 2;

    const SpotList found = SpotFinder{}.find(img, settings);
    CHECK(found.spots.size() >= 3);
    CHECK_EQ(found.header.nfast, 64);
    CHECK_EQ(found.header.nslow, 64);
    CHECK(std::fabs(found.header.pixel_mm - 0.1) < 1e-9);

    const auto tmp = std::filesystem::temp_directory_path() / "mosflm_spot_finder_test.spt";
    found.write_spt(tmp);
    const SpotList round = SpotList::read_spt(tmp);
    CHECK_EQ(round.spots.size(), found.spots.size());
    CHECK(std::fabs(round.spots[0].x_mm - found.spots[0].x_mm) < 0.05);
    std::filesystem::remove(tmp);

    // A mask (e.g. DEFPIX's beam halo) hides the peak at (16, 16).
    {
        auto mask = std::make_shared<std::vector<std::uint8_t>>(64U * 64U, 1);
        for (int r = 10; r < 23; ++r) {
            for (int c = 10; c < 23; ++c) {
                (*mask)[static_cast<std::size_t>(r) * 64U + static_cast<std::size_t>(c)] = 0;
            }
        }
        SpotFindSettings masked = settings;
        masked.mask = mask;
        const SpotList kept = SpotFinder{}.find(img, masked);
        CHECK_EQ(kept.spots.size(), found.spots.size() - 1);
        for (const auto& spot : kept.spots) {
            CHECK(std::hypot(spot.x_mm - 1.6, spot.y_mm - 1.6) > 0.5);
        }
        masked.mask = std::make_shared<std::vector<std::uint8_t>>(10U, 0);  // wrong size: ignored
        CHECK_EQ(SpotFinder{}.find(img, masked).spots.size(), found.spots.size());
    }

    // A halo (local background 200) hides a faint bump but keeps a strong spot.
    {
        auto halo = make_synthetic_frame();
        auto local = std::make_shared<std::vector<float>>(64U * 64U, -1.0f);
        for (int r = 40; r < 60; ++r) {
            for (int c = 4; c < 24; ++c) {
                const auto i = static_cast<std::size_t>(r) * 64U + static_cast<std::size_t>(c);
                halo.pixels[i] = 200;
                (*local)[i] = 200.0f;
            }
        }
        halo.pixels[50U * 64U + 8U] = 260;   // 60 over the halo: ~4σ, noise
        halo.pixels[50U * 64U + 18U] = 5000; // a real spot
        SpotFindSettings judged = settings;
        const auto near = [](const SpotList& list, double x_mm, double y_mm) {
            for (const auto& spot : list.spots) {
                if (std::hypot(spot.x_mm - x_mm, spot.y_mm - y_mm) < 0.15) {
                    return true;
                }
            }
            return false;
        };
        const SpotList plain = SpotFinder{}.find(halo, judged);
        CHECK(near(plain, 0.8, 5.0));  // the global threshold takes the bump
        judged.local_background = local;
        const SpotList local_found = SpotFinder{}.find(halo, judged);
        CHECK(!near(local_found, 0.8, 5.0));
        CHECK(near(local_found, 1.8, 5.0));
    }

    // Oracle .spt round-trip write preserves spot count.
    const std::filesystem::path insulin =
        "testdata/oracle/known_indexed/insulin_1.spt";
    if (std::filesystem::exists(insulin)) {
        const SpotList orig = SpotList::read_spt(insulin);
        const auto tmp2 =
            std::filesystem::temp_directory_path() / "mosflm_insulin_round.spt";
        orig.write_spt(tmp2);
        const SpotList again = SpotList::read_spt(tmp2);
        CHECK_EQ(again.spots.size(), orig.spots.size());
        std::filesystem::remove(tmp2);
    }

    return mosflm_test::finish("spot_finder");
}
