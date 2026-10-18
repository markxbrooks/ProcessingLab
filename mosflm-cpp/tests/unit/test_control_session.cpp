/**
 * Phase D / D5: CONTROL keyword parser unit tests.
 */

#include "mosflm/control/ControlSession.hpp"
#include "test_util.hpp"

#include <cmath>
#include <sstream>

using mosflm::control::ControlSession;

int main() {
    {
        ControlSession s;
        s.parse_line("DISTANCE 160.0");
        s.parse_line("WAVELENGTH 0.9790");
        s.parse_line("PIXEL 0.0816 0.0816");
        s.parse_line("BEAM 94.0 94.7");
        s.parse_line("CELL 78.190 78.190 78.190 90.0 90.0 90.0");
        s.parse_line("SYMMETRY 197");
        s.parse_line("MOSAIC 0.30");
        s.parse_line("MATRIX");
        s.parse_line("-0.012237 -0.000358 -0.002626");
        s.parse_line("0.002191 0.005616 -0.010974");
        s.parse_line("0.001492 -0.011185 -0.005426");
        s.parse_line("PROCESS 1 to 1 start 0.0 angle 1.0");

        const auto& st = s.state();
        CHECK(st.has_distance);
        CHECK(st.has_wavelength);
        CHECK(st.has_pixel);
        CHECK(st.has_beam);
        CHECK(st.has_cell);
        CHECK(st.has_symmetry);
        CHECK(st.has_matrix);
        CHECK(st.process_pending);
        CHECK_EQ(st.symmetry_number, 197);
        CHECK(std::fabs(st.cell.a - 78.190) < 1e-6);
        CHECK(std::fabs(st.amat(0, 0) + 0.012237) < 1e-6);
        CHECK(std::fabs(st.amat(2, 2) + 0.005426) < 1e-6);
        CHECK_EQ(st.process_first, 1);
        CHECK_EQ(st.process_last, 1);
        CHECK(std::fabs(st.process_angle_deg - 1.0) < 1e-12);
    }

    {
        ControlSession s;
        s.parse_line("AMAT -0.01 0 0  0 -0.01 0  0 0 -0.01");
        CHECK(s.state().has_matrix);
        CHECK(std::fabs(s.state().amat(1, 1) + 0.01) < 1e-12);
    }

    {
        // Unknown keywords ignored; EXIT sets flag.
        ControlSession s;
        s.parse_line("GENFILE gen.dat");
        s.parse_line("autoindex dps threshold 20");
        s.parse_line("EXIT");
        CHECK(s.exit_requested());
        CHECK(!s.state().has_matrix);
    }

    {
        // Real Fortran script parses CELL / PROCESS without throwing.
        ControlSession s;
        s.parse_file("testdata/images/known_indexed/ipmosflm.in");
        CHECK(s.state().has_cell);
        CHECK(s.state().has_symmetry);
        CHECK_EQ(s.state().symmetry_number, 197);
        CHECK(s.exit_requested());
        // No MATRIX in Fortran script — orientation comes from AUTOINDEX.
        CHECK(!s.state().has_matrix);
    }

    {
        ControlSession s;
        s.parse_file("testdata/images/known_indexed/cpp_process.in");
        CHECK(s.state().has_matrix);
        CHECK(s.state().has_cell);
        CHECK_EQ(s.state().symmetry_number, 197);
        CHECK(s.exit_requested());
        // GO was not executed by parse_file — only keyword state.
        CHECK(s.state().process_pending);
    }

    return mosflm_test::finish("control_session");
}
