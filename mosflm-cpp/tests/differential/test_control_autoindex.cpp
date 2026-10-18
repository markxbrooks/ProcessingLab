/**
 * Phase E / E1: CONTROL AUTOINDEX → CELL/MATRIX via registered DPS runner.
 */

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/indexing/ControlAutoIndex.hpp"
#include "test_util.hpp"

#include <cmath>
#include <filesystem>
#include <sstream>

using mosflm::control::ControlSession;

int main() {
    mosflm::indexing::register_control_autoindex();

    const std::filesystem::path spt =
        "testdata/oracle/known_indexed/insulin_1.spt";
    CHECK_MSG(std::filesystem::exists(spt), "missing insulin_1.spt");

    ControlSession session;
    session.parse_line("DIRECTORY testdata/oracle/known_indexed");
    session.parse_line("DISTANCE 160.0");
    session.parse_line("WAVELENGTH 0.9790");
    session.parse_line("PIXEL 0.0816 0.0816");
    session.parse_line("BEAM 94.0 94.7");
    session.parse_line("DMAX 320");
    session.parse_line(
        "AUTOINDEX DPS THRESHOLD 20 HKLDEV 0.3 NVECT 30 NVREF 30 "
        "SOLUTION 12 FILE insulin_1.spt"
    );

    std::istringstream go("GO\n");
    // Re-parse pending via go() directly:
    const auto r = session.go();
    CHECK(r.ran_autoindex);
    CHECK(session.state().has_cell);
    CHECK(session.state().has_matrix);
    CHECK(session.done_run());
    CHECK(std::fabs(session.state().cell.a - 78.19) < 2.0);

    return mosflm_test::finish("control_autoindex");
}
