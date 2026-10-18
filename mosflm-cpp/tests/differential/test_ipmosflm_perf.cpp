/**
 * Wall-time: ProcessingLab (C++) vs F77 ``ipmosflm`` on known-indexed PROCESS.
 *
 * Runs ``testdata/images/known_indexed/cpp_process.in`` (C++) and
 * ``ipmosflm_process.in`` (F77) against the CCP4 insulin SMV frame, times
 * ``MOSFLM_PERF_REPEATS`` (default 3) iterations each, prints median wall
 * times, and exits 0.
 *
 * Soft-skips (exit 0) when either binary or the image is missing so CI
 * without CCP4 / oracle still passes. Override with:
 *   MOSFLM_PROCESSINGLAB — C++ CONTROL binary
 *   MOSFLM_IPMOSFLM       — F77 ``ipmosflm`` / ``ipmosflm.oracle``
 *   MOSFLM_IMAGE          — path to ``insulin_1_001.img``
 *   MOSFLM_PERF_REPEATS   — timed iterations per side (default 3)
 */

#include "test_util.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

void skip_pass(const std::string& why) {
    std::cout << "ipmosflm_perf: SKIP (" << why << ")\n";
    std::exit(0);
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

fs::path find_integrationlab() {
    if (const char* env = std::getenv("MOSFLM_PROCESSINGLAB")) {
        return env;
    }
    if (const char* env = std::getenv("MOSFLM_THIN_IPMOSFLM")) {
        return env;
    }
    const char* rels[] = {
        "build/ProcessingLab",
        "build/ProcessingLab.exe",
        "./ProcessingLab",
        "./ProcessingLab.exe",
        "../bin/ProcessingLab",
        "../bin/ProcessingLab.exe",
        "bin/ProcessingLab",
    };
    for (const char* rel : rels) {
        if (fs::exists(rel)) {
            return rel;
        }
    }
    return {};
}

fs::path find_ipmosflm() {
    if (const char* env = std::getenv("MOSFLM_IPMOSFLM")) {
        return env;
    }
    const char* candidates[] = {
        "../bin/ipmosflm.oracle",
        "../bin/ipmosflm",
        "bin/ipmosflm.oracle",
        "bin/ipmosflm",
        "/Applications/ccp4-9/bin/ipmosflm",
        "/apps/ccp4/ccp4-9/bin/ipmosflm",
    };
    for (const char* c : candidates) {
        if (fs::exists(c)) {
            return c;
        }
    }
    if (const char* ccp4 = std::getenv("CCP4")) {
        const auto p = fs::path(ccp4) / "bin/ipmosflm";
        if (fs::exists(p)) {
            return p;
        }
    }
    return {};
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
    const char* fallbacks[] = {
        "/Applications/ccp4-9/examples/data/insulin_1_001.img",
        "/apps/ccp4/ccp4-9/examples/data/insulin_1_001.img",
    };
    for (const char* c : fallbacks) {
        if (fs::exists(c)) {
            return c;
        }
    }
    return {};
}

int repeats() {
    if (const char* env = std::getenv("MOSFLM_PERF_REPEATS")) {
        const int n = std::atoi(env);
        if (n > 0) {
            return n;
        }
    }
    return 3;
}

std::string dyld_prefix() {
#if defined(__APPLE__)
    const char* existing = std::getenv("DYLD_LIBRARY_PATH");
    std::string path =
        "/Applications/ccp4-9/lib:/opt/homebrew/opt/hdf5/lib";
    if (const char* ccp4 = std::getenv("CCP4")) {
        path = (fs::path(ccp4) / "lib").string() + ":" + path;
    }
    if (existing != nullptr && existing[0] != '\0') {
        path = std::string(existing) + ":" + path;
    }
    return "DYLD_LIBRARY_PATH=" + path + " ";
#else
    return {};
#endif
}

double median_ms(std::vector<double> samples) {
    if (samples.empty()) {
        return 0.0;
    }
    std::sort(samples.begin(), samples.end());
    const size_t n = samples.size();
    if (n % 2 == 1) {
        return samples[n / 2];
    }
    return 0.5 * (samples[n / 2 - 1] + samples[n / 2]);
}

using Clock = std::chrono::steady_clock;

double time_ms(const std::string& cmd) {
    const auto t0 = Clock::now();
    const int rc = run_cmd(cmd);
    const auto t1 = Clock::now();
    CHECK_EQ(rc, 0);
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

}  // namespace

int main() {
    const fs::path cpp_bin = find_integrationlab();
    const fs::path f77_bin = find_ipmosflm();
    const fs::path img = find_insulin_img();
    const fs::path cpp_script =
        "testdata/images/known_indexed/cpp_process.in";
    const fs::path f77_script =
        "testdata/images/known_indexed/ipmosflm_process.in";
    const fs::path mat =
        "testdata/images/known_indexed/insulin_1.mat";

    if (cpp_bin.empty() || !fs::exists(cpp_bin)) {
        skip_pass("no ProcessingLab (set MOSFLM_PROCESSINGLAB)");
    }
    if (f77_bin.empty() || !fs::exists(f77_bin)) {
        skip_pass("no F77 ipmosflm (set MOSFLM_IPMOSFLM or build fortran-oracle)");
    }
    if (img.empty() || !fs::exists(img)) {
        skip_pass("no insulin SMV (set MOSFLM_IMAGE or CCP4)");
    }
    CHECK_MSG(fs::exists(cpp_script), "missing cpp_process.in");
    CHECK_MSG(fs::exists(f77_script), "missing ipmosflm_process.in");
    CHECK_MSG(fs::exists(mat), "missing insulin_1.mat");

    const int n = repeats();
    std::cout << "ipmosflm_perf: C++=" << cpp_bin << "\n"
              << "ipmosflm_perf: F77=" << f77_bin << "\n"
              << "ipmosflm_perf: image=" << img << "\n"
              << "ipmosflm_perf: repeats=" << n << "\n";

    // --- C++ ProcessingLab PROCESS timings ---------------------------------
    std::vector<double> cpp_ms;
    cpp_ms.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const fs::path mtz =
            fs::temp_directory_path() /
            ("mosflm_perf_cpp_" + std::to_string(i) + ".mtz");
        fs::remove(mtz);
        const std::string cmd =
            cpp_bin.string() + " HKLOUT " + mtz.string() + " --image " +
            img.string() + " " + cpp_script.string() + " >/dev/null 2>/dev/null";
        const double ms = time_ms(cmd);
        CHECK_MSG(fs::exists(mtz) && fs::file_size(mtz) > 0,
                  "C++ PROCESS produced no MTZ");
        cpp_ms.push_back(ms);
        std::cout << "ipmosflm_perf: C++ iter " << (i + 1) << " " << ms
                  << " ms\n";
        fs::remove(mtz);
    }

    // --- F77 ipmosflm PROCESS timings (work dir with image + NEWMAT) -------
    std::vector<double> f77_ms;
    f77_ms.reserve(static_cast<size_t>(n));
    const std::string env_prefix = dyld_prefix();
    for (int i = 0; i < n; ++i) {
        const fs::path work =
            fs::temp_directory_path() /
            ("mosflm_perf_f77_" + std::to_string(i));
        fs::remove_all(work);
        fs::create_directories(work);
        fs::copy_file(img, work / "insulin_1_001.img");
        fs::copy_file(mat, work / "insulin_1.mat");
        fs::copy_file(f77_script, work / "ipmosflm_process.in");

        const std::string cmd =
            "cd " + work.string() + " && " + env_prefix + f77_bin.string() +
            " < ipmosflm_process.in >/dev/null 2>/dev/null";
        const double ms = time_ms(cmd);
        const fs::path out_mtz = work / "PERF.mtz";
        CHECK_MSG(fs::exists(out_mtz) && fs::file_size(out_mtz) > 0,
                  "F77 PROCESS produced no PERF.mtz");
        f77_ms.push_back(ms);
        std::cout << "ipmosflm_perf: F77 iter " << (i + 1) << " " << ms
                  << " ms\n";
        fs::remove_all(work);
    }

    const double cpp_med = median_ms(cpp_ms);
    const double f77_med = median_ms(f77_ms);
    const double speedup = (cpp_med > 0.0) ? (f77_med / cpp_med) : 0.0;

    std::cout << "ipmosflm_perf: median_C++_ms=" << cpp_med
              << " median_F77_ms=" << f77_med << " speedup_F77/C++=" << speedup
              << "\n";

    return mosflm_test::finish("ipmosflm_perf");
}
