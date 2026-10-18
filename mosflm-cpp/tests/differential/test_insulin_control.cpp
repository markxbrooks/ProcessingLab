/**
 * Phase D / D5: keyword-driven C++ PROCESS (cpp_process.in) without Fortran
 * CONTROL. Same ASU HKL recall gate as insulin_cpp_process (>= 0.85).
 */

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/crystal/MtzAdapter.hpp"
#include "mosflm/diffraction/ProcessPipeline.hpp"
#include "test_util.hpp"

#include <gemmi/symmetry.hpp>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <set>
#include <tuple>

using mosflm::control::ControlSession;
using mosflm::crystal::MtzAdapter;
using mosflm::diffraction::ProcessPipeline;

namespace {

using Hkl = std::tuple<int, int, int>;

std::filesystem::path find_insulin_img() {
    if (const char* env = std::getenv("MOSFLM_IMAGE")) {
        return env;
    }
    if (const char* ccp4 = std::getenv("CCP4")) {
        const auto p =
            std::filesystem::path(ccp4) / "examples/data/insulin_1_001.img";
        if (std::filesystem::exists(p)) {
            return p;
        }
    }
    return "/Applications/ccp4-9/examples/data/insulin_1_001.img";
}

std::set<Hkl> unique_hkl_mtz(const std::filesystem::path& mtz) {
    std::set<Hkl> out;
    for (const auto& r : MtzAdapter::read_simple(mtz)) {
        out.emplace(r.h, r.k, r.l);
    }
    return out;
}

std::set<Hkl> unique_hkl_simple(
    const std::vector<mosflm::crystal::SimpleReflection>& refs
) {
    std::set<Hkl> out;
    for (const auto& r : refs) {
        out.emplace(r.h, r.k, r.l);
    }
    return out;
}

std::set<Hkl> expand_symmetry(
    const std::set<Hkl>& in,
    const gemmi::SpaceGroup* sg
) {
    std::set<Hkl> out;
    const auto ops = sg->operations();
    for (const auto& hkl : in) {
        const gemmi::Op::Miller m{
            std::get<0>(hkl), std::get<1>(hkl), std::get<2>(hkl)
        };
        for (const auto& op : ops) {
            const gemmi::Op::Miller m2 = op.apply_to_hkl(m);
            out.emplace(m2[0], m2[1], m2[2]);
            out.emplace(-m2[0], -m2[1], -m2[2]);
        }
    }
    return out;
}

std::size_t covered_by(
    const std::set<Hkl>& query,
    const std::set<Hkl>& expanded
) {
    std::size_t n = 0;
    for (const auto& h : query) {
        if (expanded.count(h) != 0U) {
            ++n;
        }
    }
    return n;
}

}  // namespace

int main() {
    const std::filesystem::path oracle =
        "testdata/oracle/known_indexed/insulin_1_001.mtz";
    const std::filesystem::path script =
        "testdata/images/known_indexed/cpp_process.in";
    CHECK_MSG(std::filesystem::exists(oracle), "missing Fortran oracle MTZ");
    CHECK_MSG(std::filesystem::exists(script), "missing cpp_process.in");

    const auto img_path = find_insulin_img();
    CHECK_MSG(
        std::filesystem::exists(img_path),
        "missing insulin SMV (set MOSFLM_IMAGE or CCP4)"
    );

    ControlSession session;
    const auto runs = session.run_script(script, img_path);
    CHECK_MSG(!runs.empty(), "expected at least one GO");
    CHECK(runs.back().ran_process);
    const auto& result = runs.back().process;
    const auto& crystal = runs.back().crystal;

    CHECK_MSG(result.predictions.size() > 100, "too few predictions");
    CHECK_MSG(
        result.mtz_reflections.size() > 50,
        "too few integrated reflections"
    );

    const auto mtz_path =
        std::filesystem::temp_directory_path() / "mosflm_insulin_control.mtz";
    ProcessPipeline::write_mtz(mtz_path, crystal, result);

    const auto oracle_hkl = unique_hkl_mtz(oracle);
    const auto cpp_hkl = unique_hkl_simple(result.mtz_reflections);
    const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_number(197);
    CHECK_MSG(sg != nullptr, "missing I23 space group");
    const auto expanded = expand_symmetry(cpp_hkl, sg);
    const std::size_t hit = covered_by(oracle_hkl, expanded);
    const double recall =
        oracle_hkl.empty()
            ? 0.0
            : static_cast<double>(hit) / static_cast<double>(oracle_hkl.size());

    std::cout << "insulin_control: preds=" << result.predictions.size()
              << " mtz=" << result.mtz_reflections.size()
              << " oracle_hkl=" << oracle_hkl.size() << " recall=" << recall
              << '\n';

    CHECK_MSG(recall >= 0.85, "ASU HKL recall below 0.85");
    return mosflm_test::finish("insulin_control");
}
