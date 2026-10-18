#pragma once

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/geometry/UnitCell.hpp"

#include <cmath>
#include <limits>
#include <vector>

namespace mosflm::refinement {

/// One observed reflection for post-refinement (typically from MTZ).
struct PostRefObservation {
    int h{};
    int k{};
    int l{};
    /// Observed rocking angle (degrees), e.g. MOSFLM ``ROT``.
    double phi_obs_deg{};
    /// Optional weight (default 1). Ignored if ≤ 0.
    double weight{1.0};
    /// Optional observed detector pixels (NaN → φ-only for this spot).
    double xdet_obs_px{std::numeric_limits<double>::quiet_NaN()};
    double ydet_obs_px{std::numeric_limits<double>::quiet_NaN()};
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
    /// When true and observations carry XDET/YDET, include positional
    /// residuals (mm→px via ``detector`` + beam) in the cost.
    bool refine_against_detector{false};
    DetectorGeometry detector{};
    /// Direct-beam mm for ``cpp_mm_to_mosflm_xdet`` (image frame).
    double beam_x_mm{0.0};
    double beam_y_mm{0.0};
    /// Positional cost scale: ``cost += (rms_pos_px / position_px_per_deg)``.
    double position_px_per_deg{40.0};
    int max_iterations{12};
    /// Finite-difference step for missetting (degrees) and cell scale.
    double fd_misset_deg{0.05};
    double fd_cell_scale{0.001};
    /// Reject observations with |φ_pred − φ_obs| above this after a pass.
    double outlier_phi_deg{2.0};
    /// Reject positional outliers (pixels) when detector term is active.
    double outlier_pos_px{80.0};
};

struct PostRefinementResult {
    PostRefParameters parameters;
    OrientationMatrix orientation;
    UnitCell cell;
    double rms_phi_deg{0.0};
    double rms_phi_initial_deg{0.0};
    double rms_pos_px{0.0};
    double rms_pos_initial_px{0.0};
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

/**
 * Predict MOSFLM XDET/YDET pixels for one Miller index at ``phi_deg``.
 *
 * Returns false if the ray misses the detector plane.
 */
[[nodiscard]] bool predict_xdet_ydet_px(
    const OrientationMatrix& orientation,
    int h,
    int k,
    int l,
    double phi_deg,
    double wavelength_A,
    const DetectorGeometry& detector,
    double beam_x_mm,
    double beam_y_mm,
    double& xdet_px,
    double& ydet_px
);

/// Stage 2 post-refinement: cell scale, missetting, mosaic from partials.
///
/// Ports the POSTREFL parameter set (missets + cell + η) at a least-squares
/// φ-residual level, optionally with XDET/YDET positional residuals.
/// Multi-image NADD / IDXREF-style joint refinement remains for later.
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
