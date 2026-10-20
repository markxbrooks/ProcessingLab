#pragma once

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"

#include <limits>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace mosflm::diffraction {

/// One predicted spot (I/O contract from docs/FIRST_ROUTINE.md).
struct PredictedReflection {
    int h{};
    int k{};
    int l{};
    double x_mm{};
    double y_mm{};
    double phi_begin_deg{};
    double phi_end_deg{};
    bool full{true};  ///< false ⇒ partial
    /// Rocking-curve mode only: exact Ewald-crossing φ (degrees) and the
    /// fraction of the reflection's rocking curve recorded on this frame.
    /// ``phi_begin_deg``/``phi_end_deg`` then span the whole rocking curve
    /// (±3σ), not just this frame.
    double phi_calc_deg{std::numeric_limits<double>::quiet_NaN()};
    double fraction{1.0};
};

/// Acquisition rotation range (physical oscillation / still).
///
/// ``phi_end_deg == phi_start_deg`` ⇒ still.
struct ScanRange {
    double phi_start_deg{0.0};
    double phi_end_deg{0.0};
};

/// Algorithm controls for reflection prediction.
struct PredictionSettings {
    double dmin_A{1.5};
    /// Low-resolution cutoff (Å); spots with ``d > dmax`` are skipped.
    double dmax_A{200.0};
    /// Accept reflections with ``||s1| - 1/λ| ≤`` this (Å⁻¹).
    double ewald_tolerance_Ainv{0.0025};
    /// Predictor angular search increment (not physical oscillation width).
    double phi_step_deg{0.05};
    /// Crystal mosaic spread (degrees). Widens the Ewald gate by approximately
    /// ``mosaic_rad / λ`` (Å⁻¹), standing in for REEK η until a full rocking
    /// model exists. Zero ⇒ use ``ewald_tolerance_Ainv`` alone.
    double mosaic_spread_deg{0.0};
    /// Beam divergence full-angles (degrees), MOSFLM ``DIVH`` / ``DIVV``.
    /// Combined with mosaic as ``(η + max(DIVH,DIVV)) / λ`` in Å⁻¹ (REEK E0
    /// isotropic stand-in; not yet a full anisotropic rocking curve).
    double divh_deg{0.0};
    double divv_deg{0.0};
    /// Skip predictions whose d-spacing falls in MOSFLM ice-ring shells
    /// (REEK / DPSCHOOSE ``RINGTEST=1`` bands).
    bool exclude_ice_rings{false};
    /// If true (default), skip Miller indices systematically absent for
    /// ``PredictionRequest::crystal.spacegroup``.
    bool apply_systematic_absences{true};
    /// Reject spots whose observed φ interval exceeds this width (NWMAX-like).
    /// Zero disables. F77 insulin rejects ~12 "too wide in φ".
    double max_phi_width_deg{0.0};
    /// Reject near-spindle (cusp) reflections: ``|q_x| / |q| < cusp_fraction``.
    bool exclude_cusp{false};
    double cusp_fraction{0.05};
    /// Shrink the on-detector acceptance by this many pixels (REEK margin).
    bool apply_reek_detector_margin{false};
    double detector_margin_px{8.0};
    /// If false, Ewald acceptance uses ``ewald_tolerance_Ainv`` only; mosaic
    /// still affects partiality via the observed φ interval. REEK PROCESS
    /// paths set this false to avoid over-predicting.
    bool widen_ewald_with_mosaic{true};
    /// Scale applied to mosaic+divergence when widening the Ewald gate
    /// (1 = full; 0.5 ≈ REEK-like PROCESS for insulin).
    double mosaic_ewald_scale{1.0};
    /// Rocking-curve prediction (oscillation scans only): solve the Ewald
    /// crossing φc analytically and model each reflection's rocking curve as
    /// a Gaussian in φ with σ = (mosaic + max(DIVH,DIVV)) / (2.355 ζ), where
    /// ζ = |m·(ŝ1 × ŝ0)| is the inverse rotation Lorentz factor. A reflection
    /// is predicted on every frame within ``±(4σ + 2·osc)`` of φc; slices
    /// beyond the ±3σ core carry ``fraction`` < 0.001.
    bool rocking_model{false};
};

/// Cohesive inputs for stills or oscillation PREDICT.
struct PredictionRequest {
    crystal::CrystalMetadata crystal;
    OrientationMatrix orientation;
    DetectorGeometry detector;
    ScanRange scan;
    PredictionSettings settings;
};

/// Output of :meth:`ReflectionPredictor::predict`.
struct PredictionResult {
    std::vector<PredictedReflection> reflections;
};

/// Reflection prediction for stills and fine-sliced oscillations.
///
/// Frame: incident ``s0 = (0,0,1/λ)``, spindle rotation about +X, flat detector
/// at ``z = distance`` via :meth:`DetectorGeometry::lab_to_detector_mm`.
///
/// Oscillation: sample φ in ``[scan.phi_start, scan.phi_end]`` using
/// ``settings.phi_step_deg``; a spot is recorded if it meets the Ewald test on
/// the detector for any sample. ``full`` is true when the on-sphere φ interval
/// lies strictly inside the oscillation (does not touch either endpoint);
/// stills are always ``full``.
class ReflectionPredictor {
public:
    [[nodiscard]] PredictionResult predict(
        const PredictionRequest& request
    ) const;

    static void write_dump(
        const std::filesystem::path& path,
        const std::vector<PredictedReflection>& spots
    );

    [[nodiscard]] static std::vector<PredictedReflection> read_dump(
        const std::filesystem::path& path
    );

    /// Frozen analytical stills request (``predict_still.dump``).
    [[nodiscard]] static PredictionRequest analytical_still_fixture();

    /// Frozen analytical oscillation request (``predict_osc.dump``):
    /// same cell/detector as stills, ``φ ∈ [0, 0.5]°``, ``dmin = 5 Å``.
    [[nodiscard]] static PredictionRequest analytical_oscillation_fixture();
};

}  // namespace mosflm::diffraction
