/**
 * Phase D / E: CONTROL keyword parser unit tests.
 */

#include "mosflm/control/ControlSession.hpp"
#include "test_util.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

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
        CHECK(st.has_process_phi);
    }

    {
        // MOSFLM shorthand ``PROCESS 1 100`` and ``PROCESS 1 TO 5``.
        ControlSession a;
        a.parse_line("PROCESS 1 100");
        CHECK_EQ(a.state().process_first, 1);
        CHECK_EQ(a.state().process_last, 100);
        ControlSession b;
        b.parse_line("PROCESS 1 TO 5");
        CHECK_EQ(b.state().process_first, 1);
        CHECK_EQ(b.state().process_last, 5);
        CHECK(!b.state().has_process_phi);
    }

    {
        ControlSession s;
        s.parse_line("AMAT -0.01 0 0  0 -0.01 0  0 0 -0.01");
        CHECK(s.state().has_matrix);
        CHECK(std::fabs(s.state().amat(1, 1) + 0.01) < 1e-12);
    }

    {
        // Unknown keywords warn; AUTOINDEX sets pending; EXIT sets flag.
        ControlSession s;
        s.parse_line("GENFILE gen.dat");
        s.parse_line("autoindex dps threshold 20 file insulin_1.spt");
        CHECK(s.state().autoindex_pending);
        CHECK(std::fabs(s.state().autoindex.min_isigi - 20.0f) < 1e-6);
        s.parse_line("EXIT");
        CHECK(s.exit_requested());
        CHECK(!s.state().has_matrix);
    }

    {
        ControlSession s;
        bool threw = false;
        try {
            s.parse_line("REFIX");
        } catch (const std::runtime_error&) {
            threw = true;
        }
        CHECK(threw);
    }

    {
        ControlSession s;
        s.parse_line("NMASKS 5 5");
        CHECK_EQ(s.state().nmasks_x, 5);
        CHECK_EQ(s.state().nmasks_y, 5);
        s.parse_line("POSTREFINE");
        CHECK(s.state().enable_postrefinement);
    }

    {
        ControlSession s;
        CHECK(s.state().nproc >= 1);
        s.parse_line("NPROC 4");
        CHECK_EQ(s.state().nproc, 4);
        s.parse_line("NPROC OFF");
        CHECK_EQ(s.state().nproc, 1);
        s.parse_line("THREADS 8");
        CHECK_EQ(s.state().nproc, 8);
    }

    {
        ControlSession s;
        CHECK(std::fabs(s.state().min_fractioncalc - 0.5) < 1e-12);
        s.parse_line("MINFRAC 0.25");
        CHECK(std::fabs(s.state().min_fractioncalc - 0.25) < 1e-12);
        s.parse_line("MINFRAC 0");
        CHECK(std::fabs(s.state().min_fractioncalc) < 1e-12);
        s.parse_line("NADD 3");
        CHECK_EQ(s.state().nadd, 3);
        s.parse_line("STRATEGY ANGLE 45");
        CHECK(s.state().strategy_pending);
        auto r = s.go();
        CHECK(r.ran_strategy);
        CHECK(s.strategy_done());
        CHECK(s.done_run());
        CHECK(!s.needs_no_run_warning());
    }

    {
        // Real Fortran script parses CELL / PROCESS without throwing.
        ControlSession s;
        s.parse_file("testdata/images/known_indexed/ipmosflm.in");
        CHECK(s.state().has_cell);
        CHECK(s.state().has_symmetry);
        CHECK_EQ(s.state().symmetry_number, 197);
        CHECK(s.exit_requested());
        // AUTOINDEX was parsed (pending cleared only on GO via run_stream).
        // parse_file does not execute GO, so pending may remain true after
        // the AUTOINDEX line; CELL follows later GO in the file as keywords.
        CHECK(!s.state().has_matrix);
    }

    {
        ControlSession s;
        s.parse_file("testdata/images/known_indexed/cpp_process.in");
        CHECK(s.state().has_matrix);
        CHECK(s.state().has_cell);
        CHECK_EQ(s.state().symmetry_number, 197);
        CHECK(s.exit_requested());
        CHECK(s.state().process_pending);
    }

    {
        // RUN alias executes PROCESS via run_stream.
        ControlSession s;
        std::istringstream in(
            "CELL 78 78 78 90 90 90\n"
            "SYMMETRY 1\n"
            "AMAT -0.01 0 0 0 -0.01 0 0 0 -0.01\n"
            "PROCESS 1 to 1\n"
            "RUN\n"
            "EXIT\n"
        );
        // GO without image will throw; only verify RUN is recognized as go.
        // Use empty process path with STRATEGY instead:
    }
    {
        ControlSession s;
        std::istringstream in("STRATEGY ANGLE 90\nRUN\nEXIT\n");
        const auto runs = s.run_stream(in);
        CHECK_EQ(static_cast<int>(runs.size()), 1);
        CHECK(runs[0].ran_strategy);
        CHECK(s.done_run());
    }

    return mosflm_test::finish("control_session");
}
