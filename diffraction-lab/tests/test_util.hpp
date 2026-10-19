#pragma once

/// Minimal assertion helpers for standalone CTest executables (no framework).
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace mosflm_test {

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
