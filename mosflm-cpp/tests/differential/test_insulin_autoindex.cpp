/**
 * Stage 1 AUTOINDEX E0: insulin_1.spt → DPS index_fortran_ → cell / AMAT.
 *
 * Gates vs geometry.txt solution 12 (I23, cell ~78.19³).
 */

#include "mosflm/indexing/AutoIndexer.hpp"
#include "mosflm/indexing/SpotList.hpp"
#include "test_util.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>

using mosflm::indexing::AutoIndexGeometry;
using mosflm::indexing::AutoIndexSettings;
using mosflm::indexing::AutoIndexer;
using mosflm::indexing::SpotList;

int main() {
    const std::filesystem::path spt =
        "testdata/oracle/known_indexed/insulin_1.spt";
    CHECK_MSG(std::filesystem::exists(spt), "missing insulin_1.spt");

    const SpotList spots = SpotList::read_spt(spt);
    AutoIndexGeometry geo;
    geo.distance_mm = 160.0;
    geo.wavelength_A = 0.9790;
    geo.beam_x_mm = 94.0;
    geo.beam_y_mm = 94.7;
    geo.pixel_mm = 0.0816;
    // Oracle log: Maximum expected cell edge 319.93, 438 spots at I/σ≥20.
    geo.dmax_A = 320.0;

    AutoIndexSettings settings;
    settings.nvect = 30;
    settings.nvref = 30;
    settings.hkldev = 0.3f;
    settings.min_isigi = 20.0f;
    settings.prefer_system = "cI";

    const auto result = AutoIndexer{}.index(spots, geo, settings);
    std::cout << "insulin_autoindex: ok=" << result.ok
              << " ierr=" << result.ierr
              << " n_spots=" << result.n_spots_used
              << " cell=" << result.cell.describe()
              << " sys=" << result.selected_system
              << " char=" << result.selected_character << '\n';

    CHECK_MSG(result.ok, "autoindex failed");
    CHECK_EQ(result.ierr, 0);
    CHECK(result.n_spots_used > 100);

    // Cell close to regularized I23 insulin (solution 12).
    CHECK(std::fabs(result.cell.a - 78.19) < 2.0);
    CHECK(std::fabs(result.cell.b - 78.19) < 2.0);
    CHECK(std::fabs(result.cell.c - 78.19) < 2.0);
    CHECK(std::fabs(result.cell.alpha_deg - 90.0) < 3.0);
    CHECK(std::fabs(result.cell.beta_deg - 90.0) < 3.0);
    CHECK(std::fabs(result.cell.gamma_deg - 90.0) < 3.0);

    // Prefer cI when Bravais finds it; otherwise accept cubic-ish system.
    if (!result.selected_system.empty()) {
        const bool cubicish =
            result.selected_system == "cI" ||
            result.selected_system == "cP" ||
            result.selected_system == "cF" ||
            result.selected_system[0] == 'c';
        CHECK_MSG(cubicish, "expected cubic-family Bravais id, got " +
                                result.selected_system);
    }

    // AMAT not all zeros / identity.
    CHECK(result.amat.norm() > 1e-3);

    return mosflm_test::finish("insulin_autoindex");
}
