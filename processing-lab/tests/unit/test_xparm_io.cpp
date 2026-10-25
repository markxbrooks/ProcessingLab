#include "mosflm/control/ControlSession.hpp"
#include "mosflm/io/ImportedSweepGeometry.hpp"
#include "mosflm/io/XparmIO.hpp"
#include "test_util.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
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

    // Classic XDS GXPARM.XDS (76c1b, Pilatus 6M, tilted beam).
    const auto classic =
        std::filesystem::temp_directory_path() / "pl_test_classic_gxparm.xds";
    {
        std::ofstream out(classic);
        out << " XPARM.XDS    VERSION January 10, 2014\n"
            << "     1        0.0000    0.2000  0.999941 -0.000307 -0.010878\n"
            << "       0.979490       0.010032       0.003948       1.020883\n"
            << "     3     62.0381    186.9289     72.3589  90.000 100.185  90.000\n"
            << "      56.088875      18.588131      18.900845\n"
            << "     -14.252259    -109.993874     150.468033\n"
            << "      18.376116     -57.317955     -40.159519\n"
            << "         1      2463      2527    0.172000    0.172000\n"
            << "    1222.115234    1201.098267     384.959808\n"
            << "       1.000000       0.000000       0.000000\n"
            << "       0.000000       1.000000       0.000000\n"
            << "       0.000000       0.000000       1.000000\n"
            << "         1         1      2463         1      2527\n"
            << "    0.00    0.00    0.00  1.00000  0.00000  0.00000  0.00000  "
               "1.00000  0.00000\n";
    }
    const XparmGeometry c = read_xparm(classic);
    // λ comes from the 4-token row, not the 6-token row starting "1".
    CHECK_MSG(std::fabs(c.wavelength_A - 0.97949) < 1e-6, "classic wavelength");
    CHECK_EQ(c.spacegroup_number, 3);
    CHECK(c.has_amat);
    // Direct beam = ORG + F·tan(beam tilt)/q: XDS INTEGRATE.LP reports
    // 1244.12 1209.77 for this sweep.
    CHECK(c.has_direct_beam);
    CHECK_MSG(std::fabs(c.beam_fast_px - 1244.1) < 0.5, "direct beam fast");
    CHECK_MSG(std::fabs(c.beam_slow_px - 1209.8) < 0.5, "direct beam slow");
    // The beam is ~10 mrad off the detector normal: in ProcessingLab's frame
    // (beam +z) that is a detector TILT / TWIST. Refining them from zero on
    // 76c1b converged on 4.4 / −10.0 mrad; with these, predictions agree
    // with XDS INTEGRATE.HKL to 0.3 px.
    CHECK(c.has_tilt_twist);
    CHECK_MSG(std::fabs(c.tilt_mrad - 3.87) < 0.05, "classic tilt " + std::to_string(c.tilt_mrad));
    CHECK_MSG(std::fabs(c.twist_mrad + 9.83) < 0.05, "classic twist " + std::to_string(c.twist_mrad));

    mosflm::control::ControlSession s2;
    apply_xparm(s2, c);
    CHECK_MSG(
        std::fabs(s2.state().beam_x_mm - c.beam_fast_px * 0.172) < 1e-6,
        "classic BEAM x = direct beam fast"
    );
    CHECK_MSG(
        std::fabs(s2.state().beam_y_mm - c.beam_slow_px * 0.172) < 1e-6,
        "classic BEAM y = direct beam slow"
    );
    CHECK_CLOSE(s2.state().tilt_mrad, c.tilt_mrad, 1e-6);
    CHECK_CLOSE(s2.state().twist_mrad, c.twist_mrad, 1e-6);
    CHECK(!s2.state().spt_invertx);
    // XDS INTEGRATE style: orientation only, periodic refinement kept on.
    CHECK(s2.state().refine_orientation && !s2.state().refine_cell);
    CHECK(!s2.state().refine_distance && !s2.state().refine_beam);
    CHECK(s2.state().postrefine_period > 0);
    // XDS-frame AMAT has ProcessingLab's handedness (det < 0, as autoindex).
    CHECK(s2.state().amat.determinant() < 0.0);

    // ProcessingLab's own GXPARM writes NX NY as 0 0: the detector line must
    // still be read, or distance, pixel and beam fall back to defaults.
    {
        XparmGeometry own = g;
        own.nfast = 0;
        own.nslow = 0;
        own.distance_mm = 160.0;
        own.qx_mm = own.qy_mm = 0.0816;
        own.orgx_px = 1153.92;
        own.orgy_px = 1157.11;
        const auto path = std::filesystem::temp_directory_path() / "pl_test_own_gxparm.xds";
        write_xparm(path, own);
        const XparmGeometry back = read_xparm(path);
        CHECK_MSG(std::fabs(back.distance_mm - 160.0) < 1e-6, "own GXPARM distance");
        CHECK_MSG(std::fabs(back.qx_mm - 0.0816) < 1e-9, "own GXPARM pixel");
        CHECK_MSG(std::fabs(back.orgx_px - 1153.92) < 1e-3, "own GXPARM orgx");

        // Import fields (for provenance) and the orientation CONTROL lines.
        const auto im = mosflm::io::imported_from_xparm(back);
        CHECK(im.has_cell && im.has_symmetry && im.has_distance && im.has_pixel && im.has_beam);
        CHECK_MSG(std::fabs(im.beam_x_mm - 1153.92 * 0.0816) < 1e-3, "import beam x");
        const auto lines = mosflm::io::xparm_orientation_lines(back);
        CHECK_EQ(lines.size(), std::size_t{3});
        CHECK(lines[0] == "INVERTX OFF" && lines[1] == "REFINE ORIENTATION");
        CHECK(lines[2].rfind("MATRIX ", 0) == 0);
        XparmGeometry no_matrix = back;
        no_matrix.has_amat = false;
        CHECK(mosflm::io::xparm_orientation_lines(no_matrix).empty());
        std::error_code e;
        std::filesystem::remove(path, e);
    }

    std::error_code ec;
    std::filesystem::remove(tmp, ec);
    std::filesystem::remove(classic, ec);
    return mosflm_test::finish("xparm_io");
}
