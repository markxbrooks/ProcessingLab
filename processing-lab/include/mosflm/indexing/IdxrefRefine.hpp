#pragma once

#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/geometry/UnitCell.hpp"
#include "mosflm/indexing/LocalIndexAssign.hpp"
#include "mosflm/indexing/SpotList.hpp"

#include <Eigen/Dense>

#include <cstddef>
#include <limits>
#include <vector>

namespace mosflm::indexing {

/// Geometry for IDXREF spot-centroid refinement (Kabsch §2.8).
struct IdxrefGeometry {
    double distance_mm{160.0};
    double wavelength_A{0.9790};
    double beam_x_mm{94.0};
    double beam_y_mm{94.7};
    double pixel_mm{0.0816};
    double yscale{1.0};
    /// Detector omega from the ``.spt`` header (degrees). Predictions use
    /// ``ccom_from_spt_omega(omega_deg, invertx)`` so they land in the same
    /// image frame ``spots_to_rlp_xyz`` read the observations from.
    double omega_deg{0.0};
    /// INVERTX flag of the ``.spt`` (header field 6); together with
    /// ``omega_deg`` fixes the ``r_xyz`` ↔ MOSFLM A-matrix relation
    /// (:func:`mosflm_amat_to_xyz_reciprocal`).
    bool invertx{true};
};

struct IdxrefSettings {
    int max_iterations{10};
    double fd_misset_deg{0.05};
    double fd_cell_len{0.05};      ///< Å
    double fd_cell_ang{0.05};      ///< degrees
    double fd_beam_mm{0.02};
    double fd_distance_mm{0.05};
    double outlier_phi_deg{5.0};
    double outlier_pos_px{250.0};
    /// Cost scale: ``rms_pos / position_px_per_deg`` is added to ``rms_phi``
    /// (degrees). Per-frame spot φ (frame centre, σ ≈ 0.1–0.2° for 0.15°
    /// frames) is far noisier than XY (σ ≈ 1 px), so ≈10 px/deg balances
    /// the two; the former 40 let φ drive the misset and degraded XY
    /// (thaumatin: rms_xy 1.3 → 2.2 px, PROCESS predictions off by ~2.5 px).
    double position_px_per_deg{10.0};
    double max_misset_step_deg{0.35};
    /// Soft HKL regularizer weight (XY + φ are primary Kabsch residuals).
    double hkl_cost_weight{1.0};
    /// Reject refine if indexed fraction drops by more than this.
    double max_indexed_frac_drop{0.05};
    /// Refine detector distance when observation φ span ≥ this (degrees).
    double refine_distance_min_phi_span_deg{20.0};
    /// Relative bounds on distance (± fraction of start).
    double distance_frac_bound{0.05};
    /// INDEX_ORIGIN scan half-range (try dh,dk,dl ∈ [-n…n]).
    int index_origin_max{1};
};

struct IdxrefObservation {
    int h{0};
    int k{0};
    int l{0};
    double phi_obs_deg{0.0};
    double xdet_obs_px{0.0};
    double ydet_obs_px{0.0};
    double weight{1.0};
    /// Reciprocal lab vector (Å⁻¹) in AutoIndexer / MOSFLM A·hkl frame.
    /// Soft regularizer only; Kabsch §2.8 drives on X/Y/φ.
    double rx{std::numeric_limits<double>::quiet_NaN()};
    double ry{std::numeric_limits<double>::quiet_NaN()};
    double rz{std::numeric_limits<double>::quiet_NaN()};
};

struct IdxrefResult {
    bool ok{false};
    bool accepted{false};
    UnitCell cell{};
    Eigen::Matrix3d amat{Eigen::Matrix3d::Identity()};
    OrientationMatrix orientation{};
    double beam_x_mm{0.0};
    double beam_y_mm{0.0};
    double distance_mm{0.0};
    double distance_mm_initial{0.0};
    int index_origin_h{0};
    int index_origin_k{0};
    int index_origin_l{0};
    double rms_phi_initial{0.0};
    double rms_phi{0.0};
    double rms_pos_initial{0.0};
    double rms_pos{0.0};
    std::size_t n_used{0};
    std::size_t n_input{0};
    double indexed_frac_before{0.0};
    double indexed_frac_after{0.0};
};

/**
 * Kabsch INDEX_ORIGIN: try constant HKL offsets and keep the lowest cost.
 *
 * Modifies ``observations`` in place (adds the winning offset to each HKL).
 * \return Winning ``(dh,dk,dl)`` (may be 0,0,0).
 */
[[nodiscard]] Eigen::Vector3i select_index_origin(
    const UnitCell& cell0,
    const Eigen::Matrix3d& amat0_mosflm,
    std::vector<IdxrefObservation>& observations,
    const IdxrefGeometry& geometry,
    const IdxrefSettings& settings = {}
);

/**
 * XDS IDXREF §2.8: refine misset, full cell metric, beam, and (when φ span
 * is large) detector distance against ``.spt`` centroids.
 *
 * Reuses ``predict_xdet_ydet_px`` / Ewald φ from PostRefinement.
 */
[[nodiscard]] IdxrefResult refine_idxref(
    const UnitCell& cell0,
    const Eigen::Matrix3d& amat0_mosflm,
    const std::vector<IdxrefObservation>& observations,
    const IdxrefGeometry& geometry,
    const IdxrefSettings& settings = {}
);

/**
 * Build IDXREF observations from a spot list + local (or global) HKL.
 *
 * ``phi_obs`` = ``phi_end_deg`` (AutoIndexer RLP convention).
 * XDET/YDET via ``cpp_mm_to_mosflm_xdet_yscale`` so obs match
 * ``predict_xdet_ydet_px`` (image-frame mm → MOSFLM pixels).
 *
 * \param keep Indices into ``spots.spots`` matching reciprocal xyz rows.
 */
[[nodiscard]] std::vector<IdxrefObservation> build_idxref_observations(
    const SpotList& spots,
    const std::vector<std::size_t>& keep,
    const LocalIndexResult& local,
    const IdxrefGeometry& geometry,
    const std::vector<float>* x = nullptr,
    const std::vector<float>* y = nullptr,
    const std::vector<float>* z = nullptr
);

}  // namespace mosflm::indexing
