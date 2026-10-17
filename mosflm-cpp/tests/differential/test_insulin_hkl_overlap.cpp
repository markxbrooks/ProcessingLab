/**
 * Known-indexed HKL overlap: C++ PREDICT vs real ipmosflm insulin MTZ.
 *
 * Geometry from testdata/images/known_indexed/geometry.txt (CCP4 insulin_1_001).
 * Positional (XDET) mapping is still open; this test measures Miller-index
 * set overlap and probes a few A-matrix axis/sign conventions.
 */

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/crystal/MtzAdapter.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "test_util.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <set>
#include <tuple>
#include <vector>

using mosflm::OrientationMatrix;
using mosflm::UnitCell;
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

std::size_t intersection_size(const std::set<Hkl>& a, const std::set<Hkl>& b) {
    std::size_t n = 0;
    for (const auto& h : a) {
        if (b.count(h) != 0U) {
            ++n;
        }
    }
    return n;
}

Eigen::Matrix3d permute_axes(const Eigen::Matrix3d& A, int perm, int signs) {
    // perm: 0=XYZ, 1=XZY, 2=YXZ, 3=YZX, 4=ZXY, 5=ZYX
    static const int map[6][3] = {
        {0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    Eigen::Matrix3d P = Eigen::Matrix3d::Zero();
    for (int r = 0; r < 3; ++r) {
        const int c = map[perm][r];
        const double s = ((signs >> r) & 1) != 0 ? -1.0 : 1.0;
        P(r, c) = s;
    }
    return P * A;
}

PredictionRequest insulin_request(const Eigen::Matrix3d& A) {
    PredictionRequest req;
    req.crystal = CrystalMetadata::from_parameters(
        78.190, 78.190, 78.190, 90.0, 90.0, 90.0, "I 2 3"
    );
    req.orientation = OrientationMatrix(A);
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
    req.settings.dmin_A = 2.5;
    req.settings.ewald_tolerance_Ainv = 0.01;
    req.settings.phi_step_deg = 0.05;
    return req;
}

}  // namespace

int main() {
    const std::filesystem::path mtz =
        "testdata/oracle/known_indexed/insulin_1_001.mtz";
    CHECK_MSG(std::filesystem::exists(mtz), "missing insulin MTZ");

    const auto mtz_hkl = unique_hkl_mtz(mtz);
    CHECK_MSG(mtz_hkl.size() > 100, "MTZ HKL set too small");

    // Orientation Matrix [A] from geometry.txt (ipmosflm Crystal Parameters).
    Eigen::Matrix3d A;
    A << -0.012237, -0.000358, -0.002626, 0.002191, 0.005616, -0.010974,
        0.001492, -0.011185, -0.005426;

    ReflectionPredictor pred;
    std::size_t best_overlap = 0;
    int best_perm = -1;
    int best_signs = -1;
    std::size_t best_npred = 0;

    for (int perm = 0; perm < 6; ++perm) {
        for (int signs = 0; signs < 8; ++signs) {
            const Eigen::Matrix3d Ap = permute_axes(A, perm, signs);
            const auto result = pred.predict(insulin_request(Ap));
            const auto pred_hkl = unique_hkl_pred(result.reflections);
            const std::size_t ov = intersection_size(mtz_hkl, pred_hkl);
            if (ov > best_overlap) {
                best_overlap = ov;
                best_perm = perm;
                best_signs = signs;
                best_npred = pred_hkl.size();
            }
        }
    }

    const double recall =
        static_cast<double>(best_overlap) / static_cast<double>(mtz_hkl.size());
    std::cout << "insulin_hkl_overlap: mtz_unique=" << mtz_hkl.size()
              << " best_pred_unique=" << best_npred
              << " overlap=" << best_overlap << " recall=" << recall
              << " perm=" << best_perm << " signs=" << best_signs << '\n';

    // Soft gate while lab-frame mapping is unfinished (baseline ~7% recall
    // with perm=YZX / signs=2 on 2026-10-17). Tighten when XDET matches.
    CHECK_MSG(best_overlap > 0, "no HKL overlap for any A-matrix convention");
    CHECK_MSG(best_npred > 0, "predictor returned empty set");
    CHECK_MSG(recall > 0.05, "HKL recall below recorded baseline (0.05)");

    return mosflm_test::finish("insulin_hkl_overlap");
}
