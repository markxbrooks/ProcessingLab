/**
 * XDS-style completeness / quality table from unmerged PROCESS rows.
 */

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/diffraction/DataQualityStats.hpp"
#include "test_util.hpp"

#include <sstream>
#include <string>
#include <vector>

using mosflm::crystal::CrystalMetadata;
using mosflm::crystal::MosflmReflection;
using mosflm::diffraction::compute_data_quality;
using mosflm::diffraction::write_data_quality_table;

namespace {

MosflmReflection make_ref(int h, int k, int l, float i, float s) {
    MosflmReflection r;
    r.h = h;
    r.k = k;
    r.l = l;
    r.intensity = i;
    r.sigma = s;
    r.flag = 0;
    return r;
}

}  // namespace

int main() {
    const auto crystal = CrystalMetadata::from_parameters(
        78.19, 78.19, 78.19, 90.0, 90.0, 90.0, "I 2 3"
    );

    // I23 requires h+k+l even. Multiplicity ≥ 2 → R-factor / R-meas / CC.
    std::vector<MosflmReflection> refs = {
        make_ref(2, 0, 0, 100.0f, 5.0f),
        make_ref(2, 0, 0, 110.0f, 5.0f),
        make_ref(0, 2, 0, 105.0f, 5.0f),  // ASU mate of (2,0,0) under I23
        make_ref(4, 0, 0, 50.0f, 4.0f),
        make_ref(4, 0, 0, 40.0f, 4.0f),
        make_ref(3, 1, 0, 20.0f, 2.0f),
        make_ref(3, 1, 0, 24.0f, 2.0f),
    };

    const auto report =
        compute_data_quality(crystal, refs, 200.0, 1.45, 9, -3.0);
    CHECK(report.n_accepted_obs >= 5);
    CHECK(report.n_unique_accepted >= 2);
    CHECK(report.total.observed == report.n_accepted_obs);
    CHECK(report.total.possible > 0);
    CHECK(report.total.compared >= 4);
    CHECK(report.total.r_obs_pct > 0.0);
    CHECK(report.total.r_obs_pct < 100.0);
    CHECK(report.total.rmeas_pct > 0.0);
    CHECK(static_cast<int>(report.shells.size()) == 9);

    std::ostringstream os;
    write_data_quality_table(os, report);
    const std::string text = os.str();
    CHECK_MSG(
        text.find("COMPLETENESS AND QUALITY OF DATA SET") != std::string::npos,
        "missing completeness header"
    );
    CHECK_MSG(
        text.find("CC(1/2)") != std::string::npos, "missing CC(1/2) legend"
    );
    CHECK_MSG(text.find("R-meas") != std::string::npos, "missing R-meas");
    CHECK_MSG(
        text.find("NUMBER OF UNIQUE ACCEPTED REFLECTIONS") !=
            std::string::npos,
        "missing unique count footer"
    );
    CHECK_MSG(text.find("total") != std::string::npos, "missing total row");

    return mosflm_test::finish("data_quality_stats");
}
