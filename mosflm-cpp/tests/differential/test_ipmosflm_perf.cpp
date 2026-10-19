/**
 * Performance fixture: ProcessingLab ``IntegrationLab`` (C++ CONTROL) vs
 * F77 ``ipmosflm`` (CCP4 or local oracle) on the same known-indexed insulin
 * PROCESS workload.
 *
 * Soft-skip (exit 0) if either binary / image is missing or the binary is
 * the wrong architecture for this host. Override paths with:
 *   MOSFLM_INTEGRATIONLAB  — C++ binary (default: build/IntegrationLab)
 *   MOSFLM_IPMOSFLM        — F77 binary (default: $CCP4/bin/ipmosflm)
 *   MOSFLM_IMAGE           — SMV frame (default: CCP4 insulin_1_001.img)
 *   MOSFLM_PERF_REPEATS    — timed runs per binary (default: 3)
 *
 * Fixtures:
 *   testdata/images/known_indexed/cpp_process.in       — C++ CONTROL script
 *   testdata/images/known_indexed/ipmosflm_process.in  — F77 PROCESS-only
 *   testdata/images/known_indexed/insulin_1.mat        — NEWMAT for F77 MATRIX
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

bool looks_executable_here(const fs::path& bin) {
    if (!fs::exists(bin) || !fs::is_regular_file(bin)) {
        return false;
    }
    // Reject cross-arch artefacts (e.g. Mach-O under bin/ on Linux).
    const std::string probe =
        "head -c 4 \"" + bin.string() + "\" 2>/dev/null | od -An -t x1";
    FILE* pipe = popen(probe.c_str(), "r");
    if (pipe == nullptr) {
        return false;
    }
    char buf[64] = {};
    const char* got = fgets(buf, sizeof(buf), pipe);
    pclose(pipe);
    if (got == nullptr) {
        return false;
    }
    const std::string hex(buf);
#if defined(__APPLE__)
    // Mach-O magic: fe ed fa cf / cf fa ed fe (and fat variants).
    return hex.find("cf fa ed fe") != std::string::npos ||
           hex.find("fe ed fa cf") != std::string::npos ||
           hex.find("ca fe ba be") != std::string::npos ||
           hex.find("be ba fe ca") != std::string::npos;
#else
    // ELF magic: 7f 45 4c 46
    return hex.find("7f 45 4c 46") != std::string::npos;
#endif
}

fs::path find_integrationlab() {
    if (const char* env = std::getenv("MOSFLM_INTEGRATIONLAB")) {
        return env;
    }
    const char* rels[] = {
        "build/IntegrationLab",
        "build/IntegrationLab.exe",
        "./IntegrationLab",
        "./IntegrationLab.exe",
        "../build/IntegrationLab",
        "../bin/IntegrationLab",
        "bin/IntegrationLab",
    };
    for (const char* rel : rels) {
        if (looks_executable_here(rel)) {
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
        "/apps/ccp4/ccp4-9/bin/ipmosflm",
        "/Applications/ccp4-9/bin/ipmosflm",
        "../bin/ipmosflm.oracle",
        "../bin/ipmosflm",
        "bin/ipmosflm.oracle",
        "bin/ipmosflm",
    };
    for (const char* c : candidates) {
        if (looks_executable_here(c)) {
            return c;
        }
    }
    if (const char* ccp4 = std::getenv("CCP4")) {
        const auto p = fs::path(ccp4) / "bin/ipmosflm";
        if (looks_executable_here(p)) {
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
    const char* candidates[] = {
        "/apps/ccp4/ccp4-9/examples/data/insulin_1_001.img",
        "/Applications/ccp4-9/examples/data/insulin_1_001.img",
    };
    for (const char* c : candidates) {
        if (fs::exists(c)) {
            return c;
        }
    }
    return {};
}

int repeats() {
    if (const char* env = std::getenv("MOSFLM_PERF_REPEATS")) {
        try {
            const int n = std::stoi(env);
            return n > 0 ? n : 3;
        } catch (...) {
            return 3;
        }
    }
    return 3;
}

/// LD_LIBRARY_PATH prefix so CCP4 F77 ipmosflm finds ncurses.so.5 / CCP4 libs.
std::string f77_env_prefix() {
    const fs::path compat =
        fs::temp_directory_path() / "mosflm_ipmosflm_perf_compat";
    fs::create_directories(compat);
    const char* lib_dirs[] = {
        "/usr/lib/x86_64-linux-gnu",
        "/usr/lib64",
        "/lib/x86_64-linux-gnu",
    };
    for (const char* name : {"ncurses", "tinfo"}) {
        const fs::path link = compat / (std::string("lib") + name + ".so.5");
        if (fs::exists(link)) {
            continue;
        }
        for (const char* dir : lib_dirs) {
            const fs::path so6 =
                fs::path(dir) / (std::string("lib") + name + ".so.6");
            if (fs::exists(so6)) {
                std::error_code ec;
                fs::create_symlink(so6, link, ec);
                break;
            }
        }
    }

    std::string ld = compat.string() + ":/usr/lib:/usr/lib/x86_64-linux-gnu";
    if (const char* ccp4 = std::getenv("CCP4")) {
        ld += ":" + (fs::path(ccp4) / "lib").string();
    } else {
        ld += ":/apps/ccp4/ccp4-9/lib";
    }
    if (const char* existing = std::getenv("LD_LIBRARY_PATH")) {
        if (existing[0] != '\0') {
            ld += ":" + std::string(existing);
        }
    }
    std::string env = "LD_LIBRARY_PATH=" + ld;
    if (const char* ccp4 = std::getenv("CCP4")) {
        env += " CCP4=" + std::string(ccp4);
    } else {
        env += " CCP4=/apps/ccp4/ccp4-9";
    }
    return env + " ";
}

double median_ms(std::vector<double> samples) {
    if (samples.empty()) {
        return 0.0;
    }
    std::sort(samples.begin(), samples.end());
    const std::size_t n = samples.size();
    if (n % 2U == 1U) {
        return samples[n / 2U];
    }
    return 0.5 * (samples[n / 2U - 1U] + samples[n / 2U]);
}

struct TimedRun {
    double median_ms{};
    std::vector<double> samples_ms;
    fs::path mtz;
    bool ok{false};
};

TimedRun time_integrationlab(
    const fs::path& bin,
    const fs::path& img,
    const fs::path& script,
    int n
) {
    TimedRun out;
    out.mtz = fs::temp_directory_path() / "mosflm_perf_cpp.mtz";
    for (int i = 0; i < n; ++i) {
        std::error_code ec;
        fs::remove(out.mtz, ec);
        const std::string cmd = bin.string() + " HKLOUT " + out.mtz.string() +
                                " --image " + img.string() + " " +
                                script.string() + " >/dev/null 2>/dev/null";
        const auto t0 = std::chrono::steady_clock::now();
        const int rc = run_cmd(cmd);
        const auto t1 = std::chrono::steady_clock::now();
        if (rc != 0 || !fs::exists(out.mtz) || fs::file_size(out.mtz) == 0) {
            return out;
        }
        const double ms =
            std::chrono::duration<double, std::milli>(t1 - t0).count();
        out.samples_ms.push_back(ms);
    }
    out.median_ms = median_ms(out.samples_ms);
    out.ok = true;
    return out;
}

TimedRun time_ipmosflm(
    const fs::path& bin,
    const fs::path& img,
    const fs::path& script,
    const fs::path& matrix,
    int n
) {
    TimedRun out;
    const fs::path work =
        fs::temp_directory_path() / "mosflm_perf_f77_work";
    fs::create_directories(work);
    const fs::path link = work / "insulin_1_001.img";
    {
        std::error_code ec;
        fs::remove(link, ec);
        fs::create_symlink(img, link, ec);
        if (ec) {
            fs::copy_file(img, link, fs::copy_options::overwrite_existing);
        }
    }
    fs::copy_file(
        script,
        work / "ipmosflm_process.in",
        fs::copy_options::overwrite_existing
    );
    fs::copy_file(
        matrix,
        work / "insulin_1.mat",
        fs::copy_options::overwrite_existing
    );
    out.mtz = work / "PERF.mtz";
    const fs::path fallback_mtz = work / "insulin_1_001.mtz";
    const std::string env = f77_env_prefix();

    for (int i = 0; i < n; ++i) {
        std::error_code ec;
        fs::remove(out.mtz, ec);
        fs::remove(fallback_mtz, ec);
        const std::string cmd =
            "cd " + work.string() + " && " + env + bin.string() +
            " < ipmosflm_process.in >perf.log 2>&1";
        const auto t0 = std::chrono::steady_clock::now();
        const int rc = run_cmd(cmd);
        const auto t1 = std::chrono::steady_clock::now();
        if (!fs::exists(out.mtz) && fs::exists(fallback_mtz)) {
            fs::rename(fallback_mtz, out.mtz, ec);
        }
        if (rc != 0 || !fs::exists(out.mtz) || fs::file_size(out.mtz) == 0) {
            std::cerr << "ipmosflm_perf: F77 run failed (rc=" << rc
                      << "); see " << (work / "perf.log") << '\n';
            return out;
        }
        const double ms =
            std::chrono::duration<double, std::milli>(t1 - t0).count();
        out.samples_ms.push_back(ms);
    }
    out.median_ms = median_ms(out.samples_ms);
    out.ok = true;
    return out;
}

void print_samples(const char* label, const TimedRun& run) {
    std::cout << "ipmosflm_perf: " << label << " median_ms=" << run.median_ms
              << " samples_ms=[";
    for (std::size_t i = 0; i < run.samples_ms.size(); ++i) {
        if (i != 0U) {
            std::cout << ", ";
        }
        std::cout << run.samples_ms[i];
    }
    std::cout << "] mtz=" << run.mtz << '\n';
}

}  // namespace

int main() {
    const fs::path cpp_script =
        "testdata/images/known_indexed/cpp_process.in";
    const fs::path f77_script =
        "testdata/images/known_indexed/ipmosflm_process.in";
    const fs::path f77_matrix =
        "testdata/images/known_indexed/insulin_1.mat";
    CHECK_MSG(fs::exists(cpp_script), "missing cpp_process.in");
    CHECK_MSG(fs::exists(f77_script), "missing ipmosflm_process.in");
    CHECK_MSG(fs::exists(f77_matrix), "missing insulin_1.mat");

    const fs::path cpp_bin = find_integrationlab();
    if (cpp_bin.empty()) {
        skip_pass("no runnable IntegrationLab (set MOSFLM_INTEGRATIONLAB)");
    }
    const fs::path f77_bin = find_ipmosflm();
    if (f77_bin.empty()) {
        skip_pass(
            "no runnable F77 ipmosflm (set MOSFLM_IPMOSFLM, e.g. "
            "/apps/ccp4/ccp4-9/bin/ipmosflm)"
        );
    }
    const fs::path img = find_insulin_img();
    if (img.empty() || !fs::exists(img)) {
        skip_pass("no insulin SMV (set MOSFLM_IMAGE or CCP4)");
    }

    const int n = repeats();
    std::cout << "ipmosflm_perf: cpp=" << cpp_bin << " f77=" << f77_bin
              << " image=" << img << " repeats=" << n << '\n';

    const TimedRun cpp = time_integrationlab(cpp_bin, img, cpp_script, n);
    CHECK_MSG(cpp.ok, "IntegrationLab PROCESS failed");
    print_samples("IntegrationLab", cpp);

    const TimedRun f77 =
        time_ipmosflm(f77_bin, img, f77_script, f77_matrix, n);
    CHECK_MSG(f77.ok, "F77 ipmosflm PROCESS failed");
    print_samples("ipmosflm_f77", f77);

    const double speedup =
        f77.median_ms > 0.0 ? (f77.median_ms / cpp.median_ms) : 0.0;
    std::cout << "ipmosflm_perf: speedup_cpp_vs_f77=" << speedup
              << " ( >1 means C++ faster )\n";

    // Write a small machine-readable summary next to the C++ MTZ.
    {
        const fs::path summary =
            fs::temp_directory_path() / "mosflm_ipmosflm_perf.txt";
        std::ofstream out(summary);
        out << "integrationlab_ms=" << cpp.median_ms << '\n'
            << "ipmosflm_f77_ms=" << f77.median_ms << '\n'
            << "speedup_cpp_vs_f77=" << speedup << '\n'
            << "cpp_bin=" << cpp_bin << '\n'
            << "f77_bin=" << f77_bin << '\n'
            << "image=" << img << '\n'
            << "repeats=" << n << '\n';
        std::cout << "ipmosflm_perf: summary=" << summary << '\n';
    }

    return mosflm_test::finish("ipmosflm_perf");
}
