#pragma once

#include <algorithm>
#include <cstdlib>

#if defined(PROCESSINGLAB_HAS_OPENMP) && PROCESSINGLAB_HAS_OPENMP
#include <omp.h>
#endif

namespace mosflm::util {

/// Default worker count: ``PROCESSINGLAB_NPROC`` env, else OpenMP max, else 1.
/// Clamped to ``[1, 32]``.
inline int default_nproc() {
    if (const char* e = std::getenv("PROCESSINGLAB_NPROC")) {
        char* end = nullptr;
        const long v = std::strtol(e, &end, 10);
        if (end != e && v >= 1) {
            return static_cast<int>(std::min<long>(v, 32));
        }
    }
#if defined(PROCESSINGLAB_HAS_OPENMP) && PROCESSINGLAB_HAS_OPENMP
    const int m = omp_get_max_threads();
    return std::max(1, std::min(m, 32));
#else
    return 1;
#endif
}

inline void apply_nproc(int nproc) {
    const int n = std::max(1, std::min(nproc, 32));
#if defined(PROCESSINGLAB_HAS_OPENMP) && PROCESSINGLAB_HAS_OPENMP
    omp_set_num_threads(n);
    // Nested parallelism off — chunk barriers stay serial.
    omp_set_max_active_levels(1);
#else
    (void)n;
#endif
}

inline bool openmp_enabled() {
#if defined(PROCESSINGLAB_HAS_OPENMP) && PROCESSINGLAB_HAS_OPENMP
    return true;
#else
    return false;
#endif
}

}  // namespace mosflm::util
