#pragma once

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/geometry/UnitCell.hpp"

#include <vector>

namespace mosflm::refinement {

/// One observed reflection for post-refinement (typically from MTZ ``ROT``).
struct PostRefObservation {
    int h{};
    int k{};
    int l{};
    /// Observed rocking angle (degrees), e.g. MOSFLM ``ROT``.
    double phi_obs_deg{};
    /// Optional weight (default 1). Ignored if ≤ 0.
    double weight{1.0};
};

/// Parameters refined by :class:`PostRefiner` (POSTREFL Stage 2 subset).
struct PostRefParameters {
    /// Incremental missetting about lab X/Y/Z (degrees), applied as
    /// ``R(dφ) · UB₀``.
    double dphix_deg{0.0};
    double dphiy_deg{0.0};
    double dphiz_deg{0.0};
    /// Isotropic cell scale about the starting cell (1 = unchanged).
    double cell_scale{1.0};
    /// Mosaic spread (degrees) estimated from φ residuals / widths.
    double mosaic_spread_deg{0.30};
};

/// Controls for :meth:`PostRefiner::refine`. Off by default in PROCESS.
struct PostRefinementSettings {
    bool refine_missetting{true};
    bool refine_cell_scale{true};
    bool refine_mosaic{true};
    int max_iterations{12};
    /// Finite-difference step for missetting (degrees) and cell scale.
    double fd_misset_deg{0.05};
    double fd_cell_scale{0.001};
    /// Reject observations with |φ_pred − φ_obs| above this after a pass.
    double outlier_phi_deg{2.0};
};

struct PostRefinementResult {
    PostRefParameters parameters;
    OrientationMatrix orientation;
    UnitCell cell;
    double rms_phi_deg{0.0};
    double rms_phi_initial_deg{0.0};
    std::size_t n_used{0};
    std::size_t n_input{0};
    bool ok{false};
};

/// Analytic Ewald φ for spindle about +X (C++ predictor frame).
///
/// Solves ``q_y sinφ + q_z cosφ = −(λ/2)|q|²`` for ``q = UB·h``.
/// Returns up to two solutions in degrees; empty if none.
[[nodiscard]] std::vector<double> ewald_phi_solutions_deg(
    const Eigen::Vector3d& q0,
    double wavelength_A
);

/// Pick the Ewald φ closest to ``target_deg`` (or NaN if no solutions).
[[nodiscard]] double nearest_ewald_phi_deg(
    const Eigen::Vector3d& q0,
    double wavelength_A,
    double target_deg
);

/// Stage 2 post-refinement: cell scale, missetting, mosaic from partials.
///
/// Ports the POSTREFL parameter set (missets + cell + η) at a least-squares
/// φ-residual level. Not a line port of ``IDXREF``; multi-image NADD remains
/// for later. Prediction-time DIVH/DIVV gate widening is on
/// ``PredictionSettings`` (REEK E0).
class PostRefiner {
public:
    [[nodiscard]] PostRefinementResult refine(
        const crystal::CrystalMetadata& crystal,
        const OrientationMatrix& orientation0,
        double wavelength_A,
        double mosaic0_deg,
        const std::vector<PostRefObservation>& observations,
        const PostRefinementSettings& settings = {}
    ) const;

    /// Apply refined missetting + cell scale to build a new UB.
    [[nodiscard]] static OrientationMatrix apply_parameters(
        const UnitCell& cell0,
        const OrientationMatrix& orientation0,
        const PostRefParameters& p
    );
};

}  // namespace mosflm::refinement
