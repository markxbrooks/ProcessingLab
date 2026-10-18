/**
 * Phase D / D1: pure C++ insulin PROCESS (SMV → predict → sum → MTZ).
 *
 * No amalgamation / mosflm_predict_hook. ASU-aware HKL recall vs Fortran
 * oracle MTZ reuses the D0 numeric gate (>= 0.85).
 */

#include "mosflm/crystal/MtzAdapter.hpp"
#include "mosflm/diffraction/ProcessPipeline.hpp"
#include "mosflm/image/DiffractionImage.hpp"
#include "test_util.hpp"

#include <gemmi/symmetry.hpp>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <set>
#include <tuple>

using mosflm::crystal::MtzAdapter;
using mosflm::diffraction::ProcessPipeline;
using mosflm::image::SmvImageReader;

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
    CHECK_MSG(std::filesystem::exists(oracle), "missing Fortran oracle MTZ");

    const auto img_path = find_insulin_img();
    CHECK_MSG(
        std::filesystem::exists(img_path),
        "missing insulin SMV (set MOSFLM_IMAGE or CCP4)"
    );

    const auto image = SmvImageReader{}.read(img_path);
    const auto req = ProcessPipeline::known_indexed_insulin(image);
    const auto result = ProcessPipeline{}.run(req);

    CHECK_MSG(result.predictions.size() > 100, "too few predictions");
    CHECK_MSG(
        result.mtz_reflections.size() > 50,
        "too few integrated reflections"
    );

    const auto mtz_path =
        std::filesystem::temp_directory_path() / "mosflm_insulin_cpp_process.mtz";
    ProcessPipeline::write_mtz(
        mtz_path, req.crystal, result, req.detector.wavelength_A
    );

    const auto oracle_hkl = unique_hkl_mtz(oracle);
    const auto cpp_hkl = unique_hkl_simple(result.mtz_reflections);
    CHECK_MSG(oracle_hkl.size() > 100, "oracle HKL set too small");

    const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_name("I23");
    CHECK_MSG(sg != nullptr, "missing I23 space group");
    const auto cpp_expanded = expand_symmetry(cpp_hkl, sg);
    const std::size_t asu_hits = covered_by(oracle_hkl, cpp_expanded);
    const double asu_recall =
        static_cast<double>(asu_hits) / static_cast<double>(oracle_hkl.size());

    std::size_t exact = 0;
    for (const auto& h : oracle_hkl) {
        if (cpp_hkl.count(h) != 0U) {
            ++exact;
        }
    }
    const double exact_recall =
        static_cast<double>(exact) / static_cast<double>(oracle_hkl.size());

    std::cout << "insulin_cpp_process: preds=" << result.predictions.size()
              << " integrated=" << result.mtz_reflections.size()
              << " oracle_unique=" << oracle_hkl.size()
              << " cpp_unique=" << cpp_hkl.size()
              << " exact_recall=" << exact_recall
              << " asu_recall=" << asu_recall << '\n';

    // D0 gate reused as ASU-aware recall (exact indices may differ by SG/Friedel).
    CHECK_MSG(asu_recall >= 0.85, "ASU-aware HKL recall below D0 gate 0.85");

    std::filesystem::remove(mtz_path);
    return mosflm_test::finish("insulin_cpp_process");
}
