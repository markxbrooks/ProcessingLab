/**
 * Phase E / E1: CONTROL AUTOINDEX → CELL/MATRIX via registered DPS runner.
 */

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/indexing/ControlAutoIndex.hpp"
#include "test_util.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>

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

    // MULTIFRAME follows a single --image file to each frame's sibling
    // (it used to read image 1 for every frame), else TEMPLATE; HDF5
    // masters are kept as given.
    {
        namespace fs = std::filesystem;
        using mosflm::indexing::resolve_findspots_frame_path;
        const fs::path dir = fs::temp_directory_path() / "pl_frame_path_test";
        fs::create_directories(dir);
        for (const char* n : {"x_1_0001.cbf", "x_1_0002.cbf", "x_1_0037.cbf"}) {
            std::ofstream(dir / n) << "stub";
        }
        ControlSession s2;
        const std::optional<fs::path> img = dir / "x_1_0001.cbf";
        CHECK(resolve_findspots_frame_path(s2.state(), img, 37) == dir / "x_1_0037.cbf");
        CHECK(resolve_findspots_frame_path(s2.state(), img, 1) == dir / "x_1_0001.cbf");
        // Sibling missing and no TEMPLATE → error, not a silent repeat.
        bool threw = false;
        try {
            (void)resolve_findspots_frame_path(s2.state(), img, 3);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        CHECK(threw);
        s2.parse_line("DIRECTORY " + dir.string());
        s2.parse_line("TEMPLATE x_1_####.cbf");
        CHECK(resolve_findspots_frame_path(s2.state(), std::nullopt, 2) == dir / "x_1_0002.cbf");
        const std::optional<fs::path> h5 = dir / "sweep_1_master.h5";
        CHECK(resolve_findspots_frame_path(s2.state(), h5, 9) == *h5);
        std::error_code ec;
        fs::remove_all(dir, ec);
    }

    return mosflm_test::finish("control_autoindex");
}
