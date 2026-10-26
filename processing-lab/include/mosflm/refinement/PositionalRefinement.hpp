#pragma once

#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/geometry/UnitCell.hpp"

#include <gemmi/symmetry.hpp>

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace mosflm::refinement {

/// One strong spot with an observed (background-subtracted) centroid.
struct PositionalObservation {
    int h{};
    int k{};
    int l{};
    /// Rotation angle near which the spot was recorded (degrees); the model
    /// position is evaluated at the nearest Ewald crossing to this angle.
    double phi_deg{};
    /// Observed centroid in MOSFLM ``XDET``/``YDET`` pixels.
    double xdet_obs_px{};
    double ydet_obs_px{};
    double weight{1.0};
    /// Observed φ centroid (degrees), e.g. from summed rocking-curve
    /// slices. When finite, adds a φ residual and makes the rotation about
    /// the spindle observable. NaN → X/Y only.
    double phi_obs_deg{std::numeric_limits<double>::quiet_NaN()};
    /// Expected error of ``phi_obs_deg`` (degrees). When finite, the φ
    /// residual is ``Δφ / phi_sigma_deg`` (1σ counts as 1 px) instead of
    /// ``phi_weight_px_per_deg · Δφ``: wide low-ζ reflections have much
    /// less certain φ centroids and otherwise dominate the fit.
    double phi_sigma_deg{std::numeric_limits<double>::quiet_NaN()};
};

/// Options for :class:`PositionalRefiner`.
struct PositionalRefinementSettings {
    bool refine_orientation{true};
    bool refine_cell{true};
    bool refine_distance{true};
    bool refine_beam{true};
    /// Detector tilt/twist (MOSFLM ``ROTMAT2`` ray–plane model).
    bool refine_tilt_twist{true};
    int max_iterations{12};
    /// Minimum observations kept after outlier rejection.
    std::size_t min_observations{50};
    /// Reject residuals above ``max(outlier_floor_px, outlier_nsigma·rms)``.
    double outlier_nsigma{4.0};
    double outlier_floor_px{3.0};
    /// Bounds on the accepted change (safety caps).
    double max_rotation_deg{1.0};
    double max_distance_frac{0.03};
    double max_beam_shift_mm{1.0};
    double max_metric_frac{0.03};
    double max_tilt_twist_mrad{30.0};
    /// φ residual weight: 1° of φ error counts as this many pixels.
    double phi_weight_px_per_deg{10.0};
};

/// Result of a positional refinement.
struct PositionalRefinementResult {
    bool ok{false};        ///< Refinement ran and stayed within the caps.
    bool improved{false};  ///< ``ok`` and the rms residual decreased.
    OrientationMatrix orientation;  ///< Refined UB (Å⁻¹, predictor frame).
    UnitCell cell;                  ///< Symmetry-constrained refined cell.
    double distance_mm{};
    double beam_x_mm{};
    double beam_y_mm{};
    /// Refined detector tilt / twist (mrad; ``use_rotmat2`` model).
    double tilt_mrad{};
    double twist_mrad{};
    double rms_initial_px{};
    double rms_px{};
    /// Final rms of the detector residuals alone (px), and of the φ
    /// residuals (degrees) of observations that carry φ (NaN if none).
    double rms_xy_px{};
    double rms_phi_deg{std::numeric_limits<double>::quiet_NaN()};
    /// Mean detector residual, predicted − observed (px): a systematic
    /// offset between the prediction and the measured centroids.
    double mean_dx_px{};
    double mean_dy_px{};
    std::size_t n_input{};
    std::size_t n_used{};
    /// Why ``ok`` is false (empty when ok).
    std::string reason;
};

/// MOSFLM-style positional refinement against observed spot centroids.
///
/// Refines, by Levenberg–Marquardt on detector X/Y (and, where observed,
/// φ) residuals: crystal rotations about the two axes perpendicular to the
/// spindle, plus the rotation about the spindle when observations carry φ
/// (it only shifts φ and is invisible in X/Y), the
/// cell restricted to the metric allowed by the space group (via
/// :class:`MetricTensorProjector`), the crystal-to-detector distance, the
/// direct-beam position and the detector tilt/twist.
class PositionalRefiner {
public:
    /// Refine geometry.
    ///
    /// :param orientation0: Starting UB (Å⁻¹). Its metric is first projected
    ///     onto the symmetry-allowed subspace of ``spacegroup``.
    /// :param spacegroup: Space group for metric constraints (null → P1).
    /// :param wavelength_A: Wavelength (Å).
    /// :param detector: Detector geometry (``distance_mm`` is refined).
    /// :param beam_x_mm: Image-frame direct beam X (MOSFLM ``BEAM``).
    /// :param beam_y_mm: Image-frame direct beam Y.
    /// :param observations: Strong spots with observed centroids.
    /// :param settings: Options.
    /// :return: Refined geometry; ``improved`` tells whether to apply it.
    [[nodiscard]] PositionalRefinementResult refine(
        const OrientationMatrix& orientation0,
        const gemmi::SpaceGroup* spacegroup,
        double wavelength_A,
        const DetectorGeometry& detector,
        double beam_x_mm,
        double beam_y_mm,
        const std::vector<PositionalObservation>& observations,
        const PositionalRefinementSettings& settings = {}
    ) const;
};

}  // namespace mosflm::refinement
