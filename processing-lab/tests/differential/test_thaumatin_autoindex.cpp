/**
 * Thaumatin Eiger HDF5: FINDSPOTS → AUTOINDEX must recover the true lattice.
 *
 * Regression for the .spt header omega: the NeXus rotation axis is along the
 * detector fast axis, so SpotToRlp needs omega = 270° (was hard-coded 0 →
 * φ applied about slow → ~29% indexed on a 1/3-volume fractional cell).
 *
 * Soft-skip unless the master exists. Gates:
 *   - HDF5 metadata rotation_axis_omega_deg ≈ 270
 *   - FINDSPOTS .spt header omega ≈ 270
 *   - AUTOINDEX indexed ≥ 70%
 *   - cell volume within 5% of DIALS (57.94 57.95 150.30 → 5.047e5 Å³)
 *   - two shortest edges ≈ 57.9 Å
 *
 * Env overrides:
 *   MOSFLM_THAUMATIN_H5   — thaumatin_9_1_master.h5
 *   MOSFLM_THAUMATIN_LAST — last frame for FINDSPOTS (default 180; first is 2).
 *     Shorter sweeps (e.g. 61) still pick a 3× FFT3D supercell: score_basis
 *     favours denser reciprocal lattices — known basis-selection issue.
 */

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/image/ImageIO.hpp"
#include "mosflm/indexing/ControlAutoIndex.hpp"
#include "mosflm/indexing/SpotList.hpp"
#include "test_util.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

using mosflm::control::ControlSession;

namespace {

std::filesystem::path find_master() {
    if (const char* e = std::getenv("MOSFLM_THAUMATIN_H5")) {
        return e;
    }
    if (const char* home = std::getenv("HOME")) {
        return std::filesystem::path(home) / "Downloads" /
               "thaumatin_9_1_master.h5";
    }
    return "thaumatin_9_1_master.h5";
}

void skip(const std::string& why) {
    std::cout << "thaumatin_autoindex: SKIP (" << why << ")\n";
    std::exit(0);
}

double cell_volume(const mosflm::UnitCell& c) {
    constexpr double kD2R = 3.14159265358979323846 / 180.0;
    const double ca = std::cos(c.alpha_deg * kD2R);
    const double cb = std::cos(c.beta_deg * kD2R);
    const double cg = std::cos(c.gamma_deg * kD2R);
    return c.a * c.b * c.c *
           std::sqrt(1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg);
}

}  // namespace

int main() {
#if !defined(MOSFLM_HAS_HDF5)
    skip("built without HDF5");
#endif
    const std::filesystem::path master = find_master();
    if (!std::filesystem::exists(master)) {
        skip("no " + master.string());
    }
    int last = 180;
    if (const char* e = std::getenv("MOSFLM_THAUMATIN_LAST")) {
        last = std::max(3, std::atoi(e));
    }

    // 1. NeXus omega axis [-1,0,0] ∥ fast [-1,0,0] → MOSFLM omega 270°.
    const auto img = mosflm::image::read_diffraction_image(master, 2);
    CHECK_MSG(
        mosflm_test::close(img.metadata.rotation_axis_omega_deg, 270.0, 1.0),
        "omega=" + mosflm_test::fmt(img.metadata.rotation_axis_omega_deg)
    );

    mosflm::indexing::register_control_autoindex();
    mosflm::indexing::register_control_findspots();

    const std::filesystem::path spt =
        std::filesystem::temp_directory_path() / "thaumatin_autoindex.spt";

    ControlSession session;
    // Multi-frame FINDSPOTS loads frames under OpenMP; libhdf5 here is not
    // thread-safe (concurrent H5 reads crash), so stay single-threaded.
    session.parse_line("NPROC 1");
    session.parse_line("DISTANCE 199.78169135");
    session.parse_line("WAVELENGTH 0.97949913");
    session.parse_line("PIXEL 0.075 0.075");
    session.parse_line("BEAM 166.15771234 172.97150533");
    session.parse_line(
        "FINDSPOTS THRESHOLD 5 MAXSPOTS 20000 FIND 2 FILE " + spt.string() +
        " MULTIFRAME RANGE 2 " + std::to_string(last)
    );
    const auto fs = session.go(master);
    CHECK(fs.ran_findspots);

    // 2. Header omega is what SpotToRlp consumes.
    const auto list = mosflm::indexing::SpotList::read_spt(spt);
    CHECK_MSG(
        mosflm_test::close(list.header.omega_deg, 270.0, 1.0),
        "spt omega=" + mosflm_test::fmt(list.header.omega_deg)
    );

    session.parse_line(
        "AUTOINDEX DPS THRESHOLD 12 HKLDEV 0.3 NVECT 40 NVREF 40 SEPMIN 6 "
        "CLUSTER_RADIUS 2 FILE " + spt.string()
    );
    const auto ai = session.go(master);
    CHECK(ai.ran_autoindex);
    const auto& st = session.state();
    CHECK(st.has_cell);

    // 3. True lattice, not a fractional sub-cell.
    CHECK_MSG(
        st.autoindex_indexed_frac >= 0.70,
        "indexed=" + mosflm_test::fmt(st.autoindex_indexed_frac)
    );
    const double vol = cell_volume(st.cell);
    CHECK_MSG(
        std::fabs(vol / 504650.0 - 1.0) < 0.05,
        "volume=" + mosflm_test::fmt(vol)
    );
    double edges[3] = {st.cell.a, st.cell.b, st.cell.c};
    std::sort(edges, edges + 3);
    CHECK_MSG(
        mosflm_test::close(edges[0], 57.9, 1.5) &&
            mosflm_test::close(edges[1], 57.9, 1.5),
        "edges=" + mosflm_test::fmt(edges[0]) + "," + mosflm_test::fmt(edges[1])
    );

    std::error_code ec;
    std::filesystem::remove(spt, ec);
    return mosflm_test::finish("thaumatin_autoindex");
}
