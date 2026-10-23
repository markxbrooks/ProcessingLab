/**
 * 76c1b PROCESS quality with XPARM-seeded geometry (Wave 4 gates).
 *
 * Frames 1–10: FLAG=0 mean I/σ must be healthy (≥ 0.8).
 * Frames 1–50: CORRECT-lite completeness/unique must not be pathological
 * (≫ 0.2% / 125 unique from the bad-lattice XDS.md append).
 *
 * Soft-skip unless the local CBF sweep + Xia2 index XPARM exist.
 */

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/diffraction/DataQualityStats.hpp"
#include "test_util.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

using mosflm::control::ControlSession;
using mosflm::crystal::CrystalMetadata;
using mosflm::diffraction::compute_data_quality;

namespace {

void skip(const std::string& why) {
    std::cout << "76c1b_process_quality: SKIP (" << why << ")\n";
    std::exit(0);
}

std::filesystem::path find_76c1b_dir() {
    if (const char* e = std::getenv("MOSFLM_76C1B_DIR")) {
        return e;
    }
    return "/Users/brooks/Personal/xtal/building/20140929_76c1b_0002212/"
           "76c1b_0002212";
}

double mean_isigi_flag0(
    const std::vector<mosflm::crystal::MosflmReflection>& rows
) {
    double sum = 0.0;
    std::size_t n = 0;
    for (const auto& row : rows) {
        if (row.flag != 0 || !(row.sigma > 0.0f)) {
            continue;
        }
        if (!std::isfinite(row.intensity) || !std::isfinite(row.sigma)) {
            continue;
        }
        sum += static_cast<double>(row.intensity / row.sigma);
        ++n;
    }
    return (n > 0) ? (sum / static_cast<double>(n)) : 0.0;
}

}  // namespace

int main() {
    const auto dir = find_76c1b_dir();
    const auto img1 = dir / "76c1b_0002212_1_0001.cbf";
    const auto img50 = dir / "76c1b_0002212_1_0050.cbf";
    const auto xparm =
        dir / "xia2/DEFAULT/NATIVE/SWEEP1/index/XPARM.XDS";
    if (!std::filesystem::exists(img1) || !std::filesystem::exists(img50) ||
        !std::filesystem::exists(xparm)) {
        skip("missing 76c1b CBF / XPARM (set MOSFLM_76C1B_DIR)");
    }

    const std::filesystem::path work =
        std::filesystem::temp_directory_path() / "pl_76c1b_process_quality";
    std::filesystem::create_directories(work);

    // --- Short wedge I/σ gate (frames 1–10) ---
    {
        ControlSession session;
        session.set_output_directory(mosflm_test::output_dir("76c1b_process_quality_78"));
        session.parse_line("DIRECTORY " + dir.string());
        session.parse_line("TEMPLATE 76c1b_0002212_1_####.cbf");
        // Geometry from the XPARM alone (a hand-set BEAM was a workaround
        // for the CBF orientation bug and now misplaces every spot); the
        // rocking width comes from AUTOSIGMA (≈0.6–0.9° here, XDS ~1.0°).
        session.parse_line("XPARM " + xparm.string());
        session.parse_line("PARTIALSCALE ON");
        session.parse_line("AUTOSIGMA ON");
        session.parse_line(
            "HKLOUT " + (work / "pl_10.mtz").string()
        );
        session.parse_line("PROCESS 1 TO 10");
        const auto r = session.go();
        CHECK(r.ran_process);
        CHECK(r.process.mosflm_reflections.size() > 100u);

        const double mean_is =
            mean_isigi_flag0(r.process.mosflm_reflections);
        std::cout << std::fixed;
        std::cout << "76c1b_process_quality: frames 1-10 FLAG=0 mean I/σ="
                  << mean_is << " n=" << r.process.mosflm_reflections.size()
                  << '\n';
        // Mosflm on the same wedge was ~0.84; PL with fixed σ is ≥ that.
        CHECK_MSG(mean_is >= 0.80, "FLAG=0 mean I/σ below Mosflm baseline");
    }

    // --- Medium wedge CORRECT-lite tables (frames 1–50) ---
    {
        ControlSession session;
        session.set_output_directory(mosflm_test::output_dir("76c1b_process_quality_107"));
        session.parse_line("DIRECTORY " + dir.string());
        session.parse_line("TEMPLATE 76c1b_0002212_1_####.cbf");
        // Geometry from the XPARM alone (a hand-set BEAM was a workaround
        // for the CBF orientation bug and now misplaces every spot); the
        // rocking width comes from AUTOSIGMA (≈0.6–0.9° here, XDS ~1.0°).
        session.parse_line("XPARM " + xparm.string());
        session.parse_line("SYMMETRY 3");
        session.parse_line("PARTIALSCALE ON");
        session.parse_line("AUTOSIGMA ON");
        session.parse_line(
            "HKLOUT " + (work / "pl_50.mtz").string()
        );
        session.parse_line("PROCESS 1 TO 50");
        const auto r = session.go();
        CHECK(r.ran_process);
        CHECK(r.process.mosflm_reflections.size() > 500u);

        const auto& st = session.state();
        const auto crystal = CrystalMetadata::from_parameters(
            st.cell.a, st.cell.b, st.cell.c, st.cell.alpha_deg,
            st.cell.beta_deg, st.cell.gamma_deg, "P 1 21 1"
        );
        const auto quality = compute_data_quality(
            crystal, r.process.mosflm_reflections, 200.0, 1.5, 9, -3.0
        );
        std::cout << "76c1b_process_quality: frames 1-50 unique="
                  << quality.n_unique_accepted
                  << " completeness=" << quality.total.completeness_pct
                  << "% mean_isigma=" << quality.total.mean_isigma
                  << " R-meas=" << quality.total.rmeas_pct << "%\n";
        // Pathological wrong-lattice run: 0.2% / 125 unique.
        CHECK_MSG(
            quality.n_unique_accepted > 500,
            "CORRECT-lite unique still pathological"
        );
        CHECK_MSG(
            quality.total.completeness_pct > 2.0,
            "CORRECT-lite completeness still pathological"
        );
        CHECK(quality.total.mean_isigma > 0.5);
    }

    return mosflm_test::finish("76c1b_process_quality");
}
