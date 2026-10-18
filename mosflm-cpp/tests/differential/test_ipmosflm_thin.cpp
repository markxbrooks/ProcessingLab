/**
 * Phase D / D6: thin ipmosflm (C++ CONTROL) without Fortran amalgamation.
 * EXIT-only smoke + optional PROCESS when insulin SMV is available.
 */

#include "test_util.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

namespace {

fs::path find_ipmosflm() {
    if (const char* env = std::getenv("MOSFLM_THIN_IPMOSFLM")) {
        return env;
    }
    const char* rels[] = {
        "build/ipmosflm",
        "./ipmosflm",
        "../build/ipmosflm",
    };
    for (const char* rel : rels) {
        if (fs::exists(rel)) {
            return rel;
        }
    }
    return "build/ipmosflm";
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
    // POSIX: wait status; low byte 0 ⇒ exit code in high byte.
    if ((rc & 0xff) == 0) {
        return (rc >> 8) & 0xff;
    }
    return rc;
}

}  // namespace

int main() {
    const fs::path bin = find_ipmosflm();
    CHECK_MSG(fs::exists(bin), "missing thin ipmosflm binary in build/");

    // EXIT-only smoke (no amalgamation, no image).
    {
        const fs::path script =
            fs::temp_directory_path() / "mosflm_d6_exit.in";
        {
            std::ofstream out(script);
            out << "EXIT\n";
        }
        const std::string cmd =
            bin.string() + " " + script.string() + " >/dev/null";
        CHECK_EQ(run_cmd(cmd), 0);
    }

    // PROCESS path when SMV is present (same gate family as insulin_control).
    const auto img = find_insulin_img();
    if (fs::exists(img)) {
        const fs::path script =
            "testdata/images/known_indexed/cpp_process.in";
        CHECK_MSG(fs::exists(script), "missing cpp_process.in");
        const fs::path mtz =
            fs::temp_directory_path() / "mosflm_d6_thin.mtz";
        const std::string cmd = bin.string() + " HKLOUT " + mtz.string() +
                                " --image " + img.string() + " " +
                                script.string() + " >/dev/null";
        CHECK_EQ(run_cmd(cmd), 0);
        CHECK_MSG(fs::exists(mtz) && fs::file_size(mtz) > 0,
                  "thin ipmosflm did not write MTZ");
        std::cout << "ipmosflm_thin: PROCESS MTZ ok (" << mtz << ")\n";
    } else {
        std::cout << "ipmosflm_thin: SKIP PROCESS (no insulin SMV)\n";
    }

    return mosflm_test::finish("ipmosflm_thin");
}
