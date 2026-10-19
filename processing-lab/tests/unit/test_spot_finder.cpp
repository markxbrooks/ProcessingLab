/**
 * DiffractionLab-style SpotFinder + SpotList::write_spt round-trip.
 */

#include "mosflm/indexing/SpotFinder.hpp"
#include "mosflm/indexing/SpotList.hpp"
#include "test_util.hpp"

#include <cmath>
#include <filesystem>
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
