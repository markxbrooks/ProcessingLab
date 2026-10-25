/**
 * DIALS-style nearest-neighbour max cell: a 187 Å axis must lift the
 * estimate above it, and repeats of a reflection on adjacent frames merge.
 */

#include "mosflm/indexing/MaxCell.hpp"
#include "test_util.hpp"

#include <cmath>
#include <random>
#include <vector>

using mosflm::indexing::estimate_max_cell;
using mosflm::indexing::ice_ring_d_star_sq;

int main() {
    // 76c1b-like primitive cell 62 × 73 × 187 Å: a sparse, noisy RLP sample.
    std::mt19937 rng(5);
    std::uniform_real_distribution<double> keep(0.0, 1.0);
    std::normal_distribution<double> noise(0.0, 0.00005);
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> z;
    for (int h = -6; h <= 6; ++h) {
        for (int k = -6; k <= 6; ++k) {
            for (int l = -18; l <= 18; ++l) {
                if (keep(rng) > 0.4) {
                    continue;
                }
                x.push_back(static_cast<float>(h / 62.0 + noise(rng)));
                y.push_back(static_cast<float>(k / 73.0 + noise(rng)));
                z.push_back(static_cast<float>(l / 187.0 + noise(rng)));
            }
        }
    }
    const std::vector<int> group(x.size(), 0);
    const auto mc = estimate_max_cell(x, y, z, group);
    CHECK_MSG(mc.ok, "lattice estimate failed");
    CHECK_MSG(mc.max_cell_A > 187.0 && mc.max_cell_A < 400.0,
              "max_cell=" + mosflm_test::fmt(mc.max_cell_A));

    // DIALS options: a percentile and log10 bins bracket the same axis.
    mosflm::indexing::MaxCellSettings pct;
    pct.percentile = 0.95;
    const auto mp = estimate_max_cell(x, y, z, group, pct);
    CHECK_MSG(mp.ok && mp.max_cell_A > 187.0 && mp.max_cell_A < 400.0,
              "percentile max_cell=" + mosflm_test::fmt(mp.max_cell_A));
    mosflm::indexing::MaxCellSettings logbin;
    logbin.log_binning = true;
    const auto ml = estimate_max_cell(x, y, z, group, logbin);
    CHECK_MSG(ml.ok && ml.max_cell_A > 187.0 && ml.max_cell_A < 400.0,
              "log max_cell=" + mosflm_test::fmt(ml.max_cell_A));

    // Ice Ih rings: 3.90 (100), 3.67 (002), 3.44 (101) Å; (001) is absent.
    const auto rings = ice_ring_d_star_sq(3.0);
    CHECK_MSG(rings.size() >= 3, "rings=" + mosflm_test::fmt(static_cast<double>(rings.size())));
    const double expect_d[] = {3.897, 3.669, 3.441};
    for (int i = 0; i < 3 && i < static_cast<int>(rings.size()); ++i) {
        CHECK_MSG(std::fabs(1.0 / std::sqrt(rings[i]) - expect_d[i]) < 0.005,
                  "ring d=" + mosflm_test::fmt(1.0 / std::sqrt(rings[i])));
    }

    // Too few points: no estimate.
    const std::vector<float> few{0.01f, 0.02f, 0.03f};
    CHECK_MSG(!estimate_max_cell(few, few, few, {0, 0, 0}).ok, "3 points gave an estimate");

    // One reflection on two adjacent frames plus a distinct spot → 2 reflections.
    mosflm::indexing::SpotList spots;
    spots.header.pixel_mm = 0.172;
    spots.header.beam_x_mm = 100.0;
    spots.header.beam_y_mm = 100.0;
    spots.spots = {{120.0, 110.0, 0.0, 0.2, 100.0, 5.0},
                   {120.1, 110.0, 0.2, 0.4, 80.0, 5.0},
                   {140.0, 90.0, 0.0, 0.2, 90.0, 5.0}};
    const auto merged = estimate_max_cell(spots, {300.0, 1.0, 100.0, 100.0});
    CHECK_MSG(merged.n_reflections == 2, "n_reflections=" + mosflm_test::fmt(merged.n_reflections));
    CHECK_MSG(merged.n_overlaps == 0, "n_overlaps=" + mosflm_test::fmt(merged.n_overlaps));

    // Two spots of one frame chained through a neighbour: an overlap. A spot
    // within 6 px of their centroid on the same frame overlaps too.
    spots.spots = {{120.0, 110.0, 0.0, 0.2, 100.0, 5.0},
                   {120.4, 110.0, 0.2, 0.4, 80.0, 5.0},
                   {120.8, 110.0, 0.0, 0.2, 80.0, 5.0},
                   {121.35, 110.0, 0.0, 0.2, 80.0, 5.0},
                   {140.0, 90.0, 0.0, 0.2, 90.0, 5.0}};
    const auto ov = estimate_max_cell(spots, {300.0, 1.0, 100.0, 100.0});
    CHECK_MSG(ov.n_reflections == 3, "overlap n_reflections=" + mosflm_test::fmt(ov.n_reflections));
    CHECK_MSG(ov.n_overlaps == 2, "n_overlaps=" + mosflm_test::fmt(ov.n_overlaps));

    return mosflm_test::finish("max_cell");
}
