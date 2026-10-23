#pragma once

/// Minimal assertion helpers for standalone CTest executables (no framework).
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <system_error>

namespace mosflm_test {

/// Fresh, empty directory for one test's outputs: ``<build>/test_output/
/// <name>`` (CMake ``PL_TEST_OUTPUT_ROOT``), so tests never write into the
/// source tree, ``testdata`` or the image directories they read.
inline std::filesystem::path output_dir(const std::string& name) {
#ifdef PL_TEST_OUTPUT_ROOT
    const std::filesystem::path root = PL_TEST_OUTPUT_ROOT;
#else
    const std::filesystem::path root =
        std::filesystem::temp_directory_path() / "processinglab_tests";
#endif
    const std::filesystem::path dir = root / name;
    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
    std::filesystem::create_directories(dir);
    return dir;
}

inline int& failures() {
    static int n = 0;
    return n;
}

inline void report(
    bool ok,
    const char* expr,
    const char* file,
    int line,
    const std::string& extra = {}
) {
    if (ok) {
        return;
    }
    ++failures();
    std::cerr << file << ":" << line << ": CHECK failed: " << expr;
    if (!extra.empty()) {
        std::cerr << "  [" << extra << "]";
    }
    std::cerr << '\n';
}

inline bool close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol;
}

inline std::string fmt(double v) {
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    return std::string(buf);
}

inline int finish(const char* name) {
    if (failures() == 0) {
        std::cout << name << ": OK\n";
        return 0;
    }
    std::cerr << name << ": " << failures() << " failure(s)\n";
    return 1;
}

}  // namespace mosflm_test

#define CHECK(cond) ::mosflm_test::report((cond), #cond, __FILE__, __LINE__)

#define CHECK_MSG(cond, msg) \
    ::mosflm_test::report((cond), #cond, __FILE__, __LINE__, (msg))

#define CHECK_CLOSE(a, b, tol)                                                \
    ::mosflm_test::report(                                                    \
        ::mosflm_test::close((a), (b), (tol)),                                \
        #a " ~= " #b,                                                         \
        __FILE__,                                                             \
        __LINE__,                                                             \
        "got " + ::mosflm_test::fmt((a)) + " vs " + ::mosflm_test::fmt((b))  \
    )

#define CHECK_EQ(a, b)                                                        \
    ::mosflm_test::report(                                                    \
        ((a) == (b)),                                                         \
        #a " == " #b,                                                         \
        __FILE__,                                                             \
        __LINE__,                                                             \
        "got " + std::to_string((a)) + " vs " + std::to_string((b))          \
    )
