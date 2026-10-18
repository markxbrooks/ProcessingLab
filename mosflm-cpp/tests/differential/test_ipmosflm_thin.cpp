/**
 * Phase D / E: IntegrationLab (C++ CONTROL) without Fortran amalgamation.
 * Empty-stdin golden lines + EXIT smoke + optional PROCESS.
 *
 * Distinct from F77 ``ipmosflm`` (``make fortran-oracle``).
 */

#include "test_util.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace fs = std::filesystem;

namespace {

fs::path find_integrationlab() {
    if (const char* env = std::getenv("MOSFLM_INTEGRATIONLAB")) {
        return env;
    }
    // Back-compat with older test env name.
    if (const char* env = std::getenv("MOSFLM_THIN_IPMOSFLM")) {
        return env;
    }
    const char* rels[] = {
        "build/IntegrationLab",
        "build/IntegrationLab.exe",
        "./IntegrationLab",
        "./IntegrationLab.exe",
        "../build/IntegrationLab",
        "../build/IntegrationLab.exe",
    };
    for (const char* rel : rels) {
        if (fs::exists(rel)) {
            return rel;
        }
    }
    return "build/IntegrationLab";
}

fs::path find_insulin_img() {
    if (const char* env = std::getenv("MOSFLM_IMAGE")) {
        return env;
    }
    if (const char* ccp4 = std::getenv("CCP4")) {
        const auto p =
            fs::path(ccp4) / "examples/data/insulin_1_001.img";
        if (fs::exists(p)) {
            return p;
        }
    }
    return "/Applications/ccp4-9/examples/data/insulin_1_001.img";
}

int run_cmd(const std::string& cmd) {
    const int rc = std::system(cmd.c_str());
    if (rc == -1) {
        return 127;
    }
    if ((rc & 0xff) == 0) {
        return (rc >> 8) & 0xff;
    }
    return rc;
}

std::string slurp(const fs::path& path) {
    std::ifstream in(path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

}  // namespace

int main() {
    const fs::path bin = find_integrationlab();
    CHECK_MSG(fs::exists(bin), "missing IntegrationLab binary in build/");

    // Empty stdin → Mosflm-like banner + NO RUN CARD + END OF PROCESSING.
    {
        const fs::path out =
            fs::temp_directory_path() / "mosflm_e0_empty.out";
        const std::string cmd =
            bin.string() + " </dev/null >" + out.string() + " 2>/dev/null";
        CHECK_EQ(run_cmd(cmd), 0);
        const std::string text = slurp(out);
        CHECK_MSG(text.find("Mosflm version 7.4.0") != std::string::npos,
                  "missing Mosflm 7.4.0 banner");
        CHECK_MSG(text.find("NO \"RUN\" CARD GIVEN") != std::string::npos,
                  "missing NO RUN CARD warning");
        CHECK_MSG(text.find("END OF PROCESSING") != std::string::npos,
                  "missing END OF PROCESSING");
        std::cout << "integrationlab: empty-stdin golden ok\n";
    }

    // EXIT-only smoke (no amalgamation, no image) — still "ran" EXIT without RUN.
    {
        const fs::path script =
            fs::temp_directory_path() / "mosflm_d6_exit.in";
        {
            std::ofstream out(script);
            out << "EXIT\n";
        }
        const fs::path out =
            fs::temp_directory_path() / "mosflm_e0_exit.out";
        const std::string cmd = bin.string() + " " + script.string() + " >" +
                                out.string() + " 2>/dev/null";
        CHECK_EQ(run_cmd(cmd), 0);
        const std::string text = slurp(out);
        CHECK_MSG(text.find("NO \"RUN\" CARD GIVEN") != std::string::npos,
                  "EXIT-only should still warn no RUN");
    }

    // STRATEGY RUN clears no-RUN warning.
    {
        const fs::path script =
            fs::temp_directory_path() / "mosflm_e3_strategy.in";
        {
            std::ofstream out(script);
            out << "STRATEGY ANGLE 90\nRUN\nEXIT\n";
        }
        const fs::path out =
            fs::temp_directory_path() / "mosflm_e3_strategy.out";
        const std::string cmd = bin.string() + " " + script.string() + " >" +
                                out.string() + " 2>/dev/null";
        CHECK_EQ(run_cmd(cmd), 0);
        const std::string text = slurp(out);
        CHECK_MSG(text.find("NO \"RUN\" CARD GIVEN") == std::string::npos,
                  "STRATEGY RUN should not warn no RUN");
        CHECK_MSG(text.find("STRATEGY") != std::string::npos,
                  "missing STRATEGY output");
    }

    const auto img = find_insulin_img();
    if (fs::exists(img)) {
        const fs::path script =
            "testdata/images/known_indexed/cpp_process.in";
        CHECK_MSG(fs::exists(script), "missing cpp_process.in");
        const fs::path mtz =
            fs::temp_directory_path() / "mosflm_d6_thin.mtz";
        const std::string cmd = bin.string() + " HKLOUT " + mtz.string() +
                                " --image " + img.string() + " " +
                                script.string() + " >/dev/null 2>/dev/null";
        CHECK_EQ(run_cmd(cmd), 0);
        CHECK_MSG(fs::exists(mtz) && fs::file_size(mtz) > 0,
                  "IntegrationLab did not write MTZ");
        std::cout << "integrationlab: PROCESS MTZ ok (" << mtz << ")\n";
    } else {
        std::cout << "integrationlab: SKIP PROCESS (no insulin SMV)\n";
    }

    return mosflm_test::finish("integrationlab");
}
