#include "mosflm/control/ControlSession.hpp"
#include "mosflm/io/XparmIO.hpp"
#include "test_util.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>

int main() {
    using mosflm::UnitCell;
    using mosflm::io::XparmGeometry;
    using mosflm::io::apply_xparm;
    using mosflm::io::read_xparm;
    using mosflm::io::write_xparm;
    using mosflm::io::xparm_from_control;

    XparmGeometry g;
    g.spacegroup_number = 75;
    g.cell = UnitCell{58.0, 58.0, 150.0, 90.0, 90.0, 90.0};
    g.wavelength_A = 0.9795;
    g.distance_mm = 200.0;
    g.qx_mm = 0.075;
    g.qy_mm = 0.075;
    g.orgx_px = 2215.0;
    g.orgy_px = 2306.0;
    g.nfast = 4150;
    g.nslow = 4371;
    g.has_amat = true;
    g.amat << 0.01, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.01;

    const auto tmp =
        std::filesystem::temp_directory_path() / "pl_test_gxparm.xds";
    write_xparm(tmp, g);
    const XparmGeometry r = read_xparm(tmp);
    CHECK_EQ(r.spacegroup_number, 75);
    CHECK_MSG(std::fabs(r.cell.a - 58.0) < 1e-6, "cell a");
    CHECK_MSG(std::fabs(r.distance_mm - 200.0) < 1e-6, "distance");
    CHECK_MSG(std::fabs(r.qx_mm - 0.075) < 1e-9, "qx");
    CHECK_MSG(std::fabs(r.orgx_px - 2215.0) < 1e-3, "orgx");
    CHECK(r.has_amat);

    mosflm::control::ControlSession session;
    apply_xparm(session, r);
    CHECK(session.state().has_cell);
    CHECK(session.state().has_symmetry);
    CHECK(session.state().has_beam);
    CHECK(session.state().has_matrix);
    CHECK_EQ(session.state().symmetry_number, 75);
    CHECK_MSG(
        std::fabs(session.state().beam_x_mm - 2215.0 * 0.075) < 1e-3, "beam x"
    );

    const auto snap = xparm_from_control(session);
    CHECK_EQ(snap.spacegroup_number, 75);
    CHECK_MSG(std::fabs(snap.cell.c - 150.0) < 1e-6, "snap c");

    std::error_code ec;
    std::filesystem::remove(tmp, ec);
    return mosflm_test::finish("xparm_io");
}
