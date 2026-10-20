/**
 * OpenMP / NPROC smoke: util helpers + optional parallel for correctness.
 */

#include "mosflm/util/OpenMP.hpp"
#include "test_util.hpp"

#include <numeric>
#include <vector>

int main() {
    {
        const int d = mosflm::util::default_nproc();
        CHECK(d >= 1);
        CHECK(d <= 32);
        mosflm::util::apply_nproc(2);
        mosflm::util::apply_nproc(1);
    }

#if defined(PROCESSINGLAB_HAS_OPENMP) && PROCESSINGLAB_HAS_OPENMP
    {
        CHECK(mosflm::util::openmp_enabled());
        constexpr int n = 64;
        std::vector<int> out(static_cast<std::size_t>(n), 0);
        mosflm::util::apply_nproc(4);
#pragma omp parallel for num_threads(4) schedule(static)
        for (int i = 0; i < n; ++i) {
            out[static_cast<std::size_t>(i)] = i * i;
        }
        for (int i = 0; i < n; ++i) {
            CHECK_EQ(out[static_cast<std::size_t>(i)], i * i);
        }
        const int sum = std::accumulate(out.begin(), out.end(), 0);
        CHECK_EQ(sum, (n - 1) * n * (2 * n - 1) / 6);
    }
#else
    {
        CHECK(!mosflm::util::openmp_enabled());
    }
#endif

    return mosflm_test::finish("openmp_nproc");
}
