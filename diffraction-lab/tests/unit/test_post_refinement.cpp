/**
 * Post-refinement unit tests: Ewald φ solver + recovery from misset perturbation.
 */

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/geometry/Rotmat.hpp"
#include "mosflm/refinement/PostRefinement.hpp"
#include "test_util.hpp"

#include <cmath>
#include <vector>

using mosflm::OrientationMatrix;
using mosflm::UnitCell;
using mosflm::crystal::CrystalMetadata;
using mosflm::refinement::PostRefObservation;
using mosflm::refinement::PostRefiner;
using mosflm::refinement::PostRefinementSettings;
using mosflm::refinement::ewald_phi_solutions_deg;
using mosflm::refinement::nearest_ewald_phi_deg;

int main() {
    // Spindle-X Ewald: q = (0, 0.1, 0), λ=1 → rhs = −(λ/2)|q|² = −0.005
    {
        const Eigen::Vector3d q0(0.0, 0.1, 0.0);
        const auto sols = ewald_phi_solutions_deg(q0, 1.0);
        CHECK(!sols.empty());
        // qy sinφ + qz cosφ = −0.005 → 0.1 sinφ = −0.005 → sinφ = −0.05
        bool found = false;
        for (double p : sols) {
            const double s = std::sin(p * 3.14159265358979323846 / 180.0);
            if (std::fabs(s + 0.05) < 1e-6) {
                found = true;
            }
        }
        CHECK(found);
    }

    // Synthetic cubic crystal: refine back a known missetting.
    const UnitCell cell{78.19, 78.19, 78.19, 90.0, 90.0, 90.0};
    const auto crystal = CrystalMetadata::from_mosflm_cell(cell, "I 2 3");
    const OrientationMatrix true_orient =
        OrientationMatrix::from_cell_and_rotation(cell);
    const double lambda = 0.9790;

    // Build observations from true orientation (analytic φ).
    std::vector<PostRefObservation> obs;
    for (int h = -4; h <= 4; ++h) {
        for (int k = -4; k <= 4; ++k) {
            for (int l = -4; l <= 4; ++l) {
                if (h == 0 && k == 0 && l == 0) {
                    continue;
                }
                if (((h + k + l) & 1) != 0) {
                    continue;  // I-centering
                }
                const Eigen::Vector3d q0 = true_orient.miller_to_lab(h, k, l);
                const auto sols = ewald_phi_solutions_deg(q0, lambda);
                if (sols.empty()) {
                    continue;
                }
                // Prefer φ near 0.5° (insulin scan centre).
                const double phi = nearest_ewald_phi_deg(q0, lambda, 0.5);
                if (!std::isfinite(phi) || std::fabs(phi) > 5.0) {
                    continue;
                }
                obs.push_back({h, k, l, phi, 1.0});
            }
        }
    }
    CHECK(obs.size() > 20);

    // Perturb orientation by +0.4° about Y.
    const Eigen::Matrix3d Rpert = mosflm::rotmat_from_missetting(0.0, 0.4, 0.0);
    const OrientationMatrix bad(Rpert * true_orient.ub());

    PostRefinementSettings settings;
    settings.refine_missetting = true;
    settings.refine_cell_scale = true;
    settings.refine_mosaic = true;
    settings.max_iterations = 20;
    settings.outlier_phi_deg = 5.0;

    const auto result =
        PostRefiner{}.refine(crystal, bad, lambda, 0.30, obs, settings);
    CHECK(result.ok);
    CHECK(result.n_used > 10);
    CHECK(result.rms_phi_deg < result.rms_phi_initial_deg);
    CHECK(result.rms_phi_deg < 0.15);
    CHECK(std::fabs(result.parameters.dphiy_deg + 0.4) < 0.25);
    CHECK(std::fabs(result.parameters.cell_scale - 1.0) < 0.02);
    CHECK(result.parameters.mosaic_spread_deg > 0.0);
    CHECK(result.parameters.mosaic_spread_deg < 2.0);

    return mosflm_test::finish("post_refinement");
}
