/**
 * Phase D / D4: post-refinement vs insulin oracle ROT (not enabled in PROCESS).
 *
 * Starts from a perturbed Orientation Matrix [A], refines missetting / cell
 * scale / mosaic against MTZ ``ROT``, and gates residual improvement plus soft
 * agreement with the known cell / mosaic.
 */

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/geometry/Rotmat.hpp"
#include "mosflm/refinement/PostRefinement.hpp"
#include "test_util.hpp"

#include <gemmi/mtz.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

using mosflm::OrientationMatrix;
using mosflm::crystal::CrystalMetadata;
using mosflm::refinement::PostRefObservation;
using mosflm::refinement::PostRefiner;
using mosflm::refinement::PostRefinementSettings;

namespace {

Eigen::Matrix3d insulin_amat() {
    Eigen::Matrix3d A;
    A << -0.012237, -0.000358, -0.002626, 0.002191, 0.005616, -0.010974,
        0.001492, -0.011185, -0.005426;
    return A;
}

std::vector<PostRefObservation> read_oracle_rot(
    const std::filesystem::path& path
) {
    gemmi::Mtz mtz;
    mtz.read_file(path.string());
    const auto* H = mtz.column_with_label("H");
    const auto* K = mtz.column_with_label("K");
    const auto* L = mtz.column_with_label("L");
    const auto* R = mtz.column_with_label("ROT");
    const auto* F = mtz.column_with_label("FRACTIONCALC");
    CHECK_MSG(H && K && L && R, "oracle missing H/K/L/ROT");
    std::vector<PostRefObservation> out;
    out.reserve(static_cast<std::size_t>(mtz.nreflections));
    for (int i = 0; i < mtz.nreflections; ++i) {
        PostRefObservation o;
        o.h = static_cast<int>(std::lround((*H)[static_cast<std::size_t>(i)]));
        o.k = static_cast<int>(std::lround((*K)[static_cast<std::size_t>(i)]));
        o.l = static_cast<int>(std::lround((*L)[static_cast<std::size_t>(i)]));
        o.phi_obs_deg = (*R)[static_cast<std::size_t>(i)];
        // Prefer partials / near-edge (more orientation-sensitive).
        const double frac =
            F ? (*F)[static_cast<std::size_t>(i)] : 0.5;
        o.weight = (frac < 0.95) ? 1.5 : 0.5;
        if (std::fabs(o.phi_obs_deg) > 5.0) {
            continue;
        }
        out.push_back(o);
    }
    return out;
}

}  // namespace

int main() {
    const std::filesystem::path oracle =
        "testdata/oracle/known_indexed/insulin_1_001.mtz";
    CHECK_MSG(std::filesystem::exists(oracle), "missing oracle MTZ");

    const auto obs = read_oracle_rot(oracle);
    CHECK_MSG(obs.size() > 100, "too few ROT observations");

    const auto crystal = CrystalMetadata::from_parameters(
        78.190, 78.190, 78.190, 90.0, 90.0, 90.0, "I 2 3"
    );
    const OrientationMatrix good =
        OrientationMatrix::from_mosflm_orientation_matrix(insulin_amat());
    const double lambda = 0.9790;

    // Perturb like a bad indexing misset.
    const Eigen::Matrix3d Rpert =
        mosflm::rotmat_from_missetting(0.35, -0.25, 0.15);
    const OrientationMatrix bad(Rpert * good.ub());

    PostRefinementSettings settings;
    settings.refine_missetting = true;
    settings.refine_cell_scale = true;
    settings.refine_mosaic = true;
    settings.max_iterations = 25;
    settings.outlier_phi_deg = 3.0;

    const auto result =
        PostRefiner{}.refine(crystal, bad, lambda, 0.50, obs, settings);

    std::cout << "insulin_postref: n_input=" << result.n_input
              << " n_used=" << result.n_used
              << " rms0=" << result.rms_phi_initial_deg
              << " rms=" << result.rms_phi_deg
              << " dphi=(" << result.parameters.dphix_deg << ','
              << result.parameters.dphiy_deg << ','
              << result.parameters.dphiz_deg << ')'
              << " cell_scale=" << result.parameters.cell_scale
              << " mosaic=" << result.parameters.mosaic_spread_deg << '\n';

    CHECK(result.ok);
    CHECK_MSG(result.n_used > 50, "too few usable φ residuals");
    CHECK_MSG(
        result.rms_phi_deg < result.rms_phi_initial_deg,
        "φ RMS did not improve"
    );
    // Soft: recovered geometry stays near indexing solution.
    CHECK_MSG(
        std::fabs(result.parameters.cell_scale - 1.0) < 0.05,
        "cell scale drifted beyond soft gate"
    );
    CHECK_MSG(
        result.parameters.mosaic_spread_deg > 0.05 &&
            result.parameters.mosaic_spread_deg < 1.5,
        "mosaic outside soft band"
    );
    // Missetting correction should partially undo the +0.35/-0.25/+0.15 kick.
    const double misset_norm = std::hypot(
        result.parameters.dphix_deg, result.parameters.dphiy_deg,
        result.parameters.dphiz_deg
    );
    CHECK_MSG(misset_norm > 0.05, "no missetting correction applied");
    CHECK_MSG(misset_norm < 2.0, "missetting correction unphysically large");

    return mosflm_test::finish("insulin_postref");
}
