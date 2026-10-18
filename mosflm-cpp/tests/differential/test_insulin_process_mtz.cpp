/**
 * Phase C / D0: insulin PROCESS MTZ with MOSFLM_CPP_PREDICT vs Fortran oracle.
 *
 * Compares unique HKL sets. Gate: exact-index recall >= 0.85.
 */

#include "mosflm/crystal/MtzAdapter.hpp"
#include "test_util.hpp"

#include <filesystem>
#include <iostream>
#include <set>
#include <tuple>

using mosflm::crystal::MtzAdapter;

namespace {

using Hkl = std::tuple<int, int, int>;

std::set<Hkl> unique_hkl(const std::filesystem::path& mtz) {
    std::set<Hkl> out;
    for (const auto& r : MtzAdapter::read_simple(mtz)) {
        out.emplace(r.h, r.k, r.l);
    }
    return out;
}

std::size_t intersection_size(const std::set<Hkl>& a, const std::set<Hkl>& b) {
    std::size_t n = 0;
    for (const auto& h : a) {
        if (b.count(h) != 0U) {
            ++n;
        }
    }
    return n;
}

}  // namespace

int main() {
    const std::filesystem::path oracle =
        "testdata/oracle/known_indexed/insulin_1_001.mtz";
    const std::filesystem::path cpp_pred =
        "testdata/oracle/known_indexed_cpp_predict/insulin_1_001.mtz";

    CHECK_MSG(std::filesystem::exists(oracle), "missing Fortran oracle MTZ");
    CHECK_MSG(
        std::filesystem::exists(cpp_pred),
        "missing CPP-predict PROCESS MTZ (rebuild with MOSFLM_CPP_PREDICT=1)"
    );

    const auto oracle_hkl = unique_hkl(oracle);
    const auto cpp_hkl = unique_hkl(cpp_pred);
    CHECK_MSG(oracle_hkl.size() > 100, "oracle HKL set too small");
    CHECK_MSG(cpp_hkl.size() > 100, "cpp_predict HKL set too small");

    const std::size_t ov = intersection_size(oracle_hkl, cpp_hkl);
    const double recall =
        static_cast<double>(ov) / static_cast<double>(oracle_hkl.size());
    const double precision =
        static_cast<double>(ov) / static_cast<double>(cpp_hkl.size());

    std::cout << "insulin_process_mtz: oracle_unique=" << oracle_hkl.size()
              << " cpp_unique=" << cpp_hkl.size() << " overlap=" << ov
              << " recall=" << recall << " precision=" << precision << '\n';

    CHECK_MSG(ov > 0, "no HKL overlap");
    CHECK_MSG(
        recall >= 0.85,
        "exact-HKL recall below Phase C gate 0.85"
    );

    return mosflm_test::finish("insulin_process_mtz");
}
