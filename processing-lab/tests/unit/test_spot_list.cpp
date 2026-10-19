/**
 * MOSFLM .spt reader (Stage 1 AUTOINDEX scaffolding).
 */

#include "mosflm/indexing/SpotList.hpp"
#include "test_util.hpp"

#include <cmath>
#include <filesystem>

using mosflm::indexing::SpotList;

int main() {
    const std::filesystem::path path =
        "testdata/oracle/known_indexed/insulin_1.spt";
    CHECK_MSG(std::filesystem::exists(path), "missing insulin_1.spt");

    const SpotList list = SpotList::read_spt(path);
    CHECK_EQ(list.header.nfast, 2304);
    CHECK_EQ(list.header.nslow, 2304);
    CHECK(std::fabs(list.header.pixel_mm - 0.0816) < 1e-4);
    CHECK(std::fabs(list.header.beam_x_mm - 94.0) < 1e-3);
    CHECK(std::fabs(list.header.beam_y_mm - 94.7) < 1e-3);
    CHECK(list.spots.size() >= 500);
    CHECK(list.spots.size() <= 700);

    const auto tmp =
        std::filesystem::temp_directory_path() / "mosflm_spot_list_write.spt";
    list.write_spt(tmp);
    const SpotList again = SpotList::read_spt(tmp);
    CHECK_EQ(again.spots.size(), list.spots.size());
    CHECK(std::fabs(again.header.beam_x_mm - list.header.beam_x_mm) < 1e-3);
    std::filesystem::remove(tmp);

    return mosflm_test::finish("spot_list");
}
