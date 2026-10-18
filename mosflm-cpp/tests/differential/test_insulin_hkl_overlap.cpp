/**
 * Known-indexed HKL overlap: C++ PREDICT vs real ipmosflm insulin MTZ.
 *
 * Applies :func:`mosflm_amat_to_cpp_ub` to the printed Orientation Matrix [A].
 * Exact Miller indices may differ by space-group / Friedel choice; the gate is
 * ASU-aware recall (every MTZ HKL matches a predicted equivalent).
 */

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/crystal/MtzAdapter.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/geometry/LabFrame.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "test_util.hpp"

#include <gemmi/symmetry.hpp>

#include <filesystem>
#include <iostream>
#include <set>
#include <tuple>
#include <vector>

using mosflm::OrientationMatrix;
using mosflm::crystal::CrystalMetadata;
using mosflm::crystal::MtzAdapter;
using mosflm::diffraction::PredictionRequest;
using mosflm::diffraction::ReflectionPredictor;

namespace {

using Hkl = std::tuple<int, int, int>;

std::set<Hkl> unique_hkl_mtz(const std::filesystem::path& mtz) {
    std::set<Hkl> out;
    for (const auto& r : MtzAdapter::read_simple(mtz)) {
        out.emplace(r.h, r.k, r.l);
    }
    return out;
}

std::set<Hkl> unique_hkl_pred(
    const std::vector<mosflm::diffraction::PredictedReflection>& spots
) {
    std::set<Hkl> out;
    for (const auto& s : spots) {
        out.emplace(s.h, s.k, s.l);
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

std::size_t intersection_size(const std::set<Hkl>& a, const std::set<Hkl>& b) {
    std::size_t n = 0;
    for (const auto& h : a) {
        if (b.count(h) != 0U) {
            ++n;
        }
    }
    return n;
}

PredictionRequest insulin_request(const Eigen::Matrix3d& Ub) {
    PredictionRequest req;
    req.crystal = CrystalMetadata::from_parameters(
        78.190, 78.190, 78.190, 90.0, 90.0, 90.0, "I 2 3"
    );
    req.orientation = OrientationMatrix(Ub);
    req.detector.distance_mm = 160.0;
    req.detector.beam_x_mm = 0.0;
    req.detector.beam_y_mm = 0.0;
    req.detector.pixel_size_x_mm = 0.0816;
    req.detector.pixel_size_y_mm = 0.0816;
    req.detector.nfast = 2304;
    req.detector.nslow = 2304;
    req.detector.wavelength_A = 0.9790;
    req.scan.phi_start_deg = 0.0;
    req.scan.phi_end_deg = 1.0;
    req.settings.dmin_A = 1.45;
    req.settings.ewald_tolerance_Ainv = 0.01;
    req.settings.phi_step_deg = 0.05;
    return req;
}

Eigen::Matrix3d insulin_amat_mosflm() {
    Eigen::Matrix3d A;
    A << -0.012237, -0.000358, -0.002626, 0.002191, 0.005616, -0.010974,
        0.001492, -0.011185, -0.005426;
    return A;
}

}  // namespace

int main() {
    const std::filesystem::path mtz =
        "testdata/oracle/known_indexed/insulin_1_001.mtz";
    CHECK_MSG(std::filesystem::exists(mtz), "missing insulin MTZ");

    const auto mtz_hkl = unique_hkl_mtz(mtz);
    CHECK_MSG(mtz_hkl.size() > 100, "MTZ HKL set too small");

    const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_name("I 2 3");
    CHECK_MSG(sg != nullptr, "space group I 2 3 not found");

    const Eigen::Matrix3d Ub =
        OrientationMatrix::from_mosflm_orientation_matrix(insulin_amat_mosflm())
            .ub();

    ReflectionPredictor pred;
    const auto result = pred.predict(insulin_request(Ub));
    const auto pred_hkl = unique_hkl_pred(result.reflections);
    const auto pred_exp = expand_symmetry(pred_hkl, sg);

    const std::size_t exact_ov = intersection_size(mtz_hkl, pred_hkl);
    const std::size_t asu_ov = intersection_size(mtz_hkl, pred_exp);
    const double exact_recall =
        static_cast<double>(exact_ov) / static_cast<double>(mtz_hkl.size());
    const double asu_recall =
        static_cast<double>(asu_ov) / static_cast<double>(mtz_hkl.size());

    std::cout << "insulin_hkl_overlap: mtz_unique=" << mtz_hkl.size()
              << " pred_unique=" << pred_hkl.size()
              << " exact_overlap=" << exact_ov
              << " exact_recall=" << exact_recall
              << " asu_overlap=" << asu_ov << " asu_recall=" << asu_recall
              << '\n';

    CHECK_MSG(pred_hkl.size() > 0, "predictor returned empty set");
    CHECK_MSG(exact_ov > 0, "no exact HKL overlap");
    // Lab-frame map is correct when every integrated HKL is a predicted
    // symmetry / Friedel equivalent (exact indices may differ).
    CHECK_MSG(
        asu_recall > 0.95,
        "ASU-aware HKL recall below 0.95 — lab-frame map may be wrong"
    );

    return mosflm_test::finish("insulin_hkl_overlap");
}
