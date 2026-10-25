#pragma once

#include "mosflm/crystal/MtzAdapter.hpp"

#include <optional>
#include <vector>

namespace mosflm::diffraction {

/// Options for :func:`sum_partial_slices`.
struct PartialSummationSettings {
    /// Slices of the same original HKL whose calculated φ differ by more than
    /// this (degrees) belong to different reflections (e.g. 180° apart).
    double max_phic_gap_deg{1.0};
    /// Slices below this fraction do not contribute flag bits (weak tails).
    double flag_min_fraction{0.2};
    /// Slices below this fraction (beyond the ±3σ core) are left out of the
    /// summed intensity but still used for the observed φ statistics.
    double core_min_fraction{1e-3};
};

/// MOSFLM-style partials summation of rocking-curve slices.
///
/// Groups rows by original Miller index (``h0,k0,l0``) and calculated φ
/// (``rot``), then writes one row per reflection: intensities summed and
/// variances added over the core slices, ``FRACTIONCALC`` = summed fraction (capped at 1),
/// ``BATCH``/``XDET``/``YDET``/``LP`` from the slice holding the largest
/// fraction, ``phi_obs`` and the observed centroid intensity-weighted over
/// slices. Overlap / background / peak-ratio flags are OR-ed over the main
/// slices; the negative-profile flag is recomputed on the sum. Reflections
/// with no core slice (core frames masked or rejected) are dropped.
///
/// :param slices: Per-frame rows (order irrelevant).
/// :param settings: Options.
/// :return: Summed reflections sorted by batch.
[[nodiscard]] std::vector<crystal::MosflmReflection> sum_partial_slices(
    std::vector<crystal::MosflmReflection> slices,
    const PartialSummationSettings& settings = {}
);

/// Result of :func:`fit_rocking_width`.
struct RockingWidthFit {
    /// Fitted rocking width (FWHM before the Lorentz factor, degrees): the
    /// predictor's ``MOSAIC + DIV``.
    double width_deg{};
    /// Strong, complete reflections used.
    std::size_t n_reflections{};
};

/// Fit one rocking width to the per-frame intensity of strong reflections.
///
/// Each reflection's slices are compared with the Gaussian model the
/// predictor uses, ``σφ = W / (2.355 ζ)``; ζ comes from the slice ``width``
/// (``6σφ`` at ``width_now_deg``). The amplitude is free per reflection, so
/// a rocking curve truncated by the measurement reach is modelled rather
/// than mistaken for a narrow one, and background noise in the outer slices
/// enters symmetrically instead of inflating a second moment. A free
/// per-slice baseline absorbs a background bias, which would otherwise read
/// as heavier tails and widen the fit at every re-estimate. Minimises the
/// summed, per-reflection normalised misfit over ``W`` in [0.02°, 3°].
///
/// :param slices: Per-frame rocking-mode rows (``rot`` = φ calc, ``phi_obs``
///     = frame mid-point).
/// :param width_now_deg: Width the slices were predicted with.
/// :param osc_deg: Oscillation width of one frame.
/// :param settings: Grouping options (as for :func:`sum_partial_slices`).
/// :param min_reflections: Fewer usable reflections → no fit.
/// :return: The fit, or ``std::nullopt`` if too few reflections qualify.
[[nodiscard]] std::optional<RockingWidthFit> fit_rocking_width(
    const std::vector<crystal::MosflmReflection>& slices,
    double width_now_deg,
    double osc_deg,
    const PartialSummationSettings& settings = {},
    std::size_t min_reflections = 40
);

}  // namespace mosflm::diffraction
