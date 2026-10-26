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
/// slices. Background and peak-ratio flags are OR-ed over the main slices;
/// the negative-profile flag is recomputed on the sum and set only when the
/// summed profile intensity is below −5σ. Reflections with no core slice
/// (core frames masked or rejected) are dropped.
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
    /// Median reflecting-range standard deviation σ_M (degrees), XDS
    /// ``REFLECTING_RANGE_E.S.D.``. The predictor uses ``σφ = σ_M / ζ``.
    double sigma_m_deg{};
    /// Reflections whose core contributed to the median.
    std::size_t n_reflections{};
};

/// Estimate σ_M from the core of each strong reflection's rocking curve.
///
/// For each reflection, ζ is recovered from the latest slice (``width`` is
/// ``6 σφ`` predicted at ``sigma_m_now_deg``). Slices with
/// ``|ε| = |ζ (φ_obs − φ_calc)|`` above 1° are ignored, so the sample does
/// not grow when the prediction reach does. A constant background is the
/// median intensity beyond 2.5σ of a first pass (and at least 0.75° from
/// φ calc). The reflecting-range variance is
/// the intensity-weighted second moment of ε minus the frame-bin term
/// ``(ζ Δφ)² / 12``. σ_M is the median of those standard deviations over
/// reflections with at least three frames in the cap and summed I/σ ≥ 10.
/// Spatially overlapped slices are kept: their φ extent is still a width.
///
/// :param slices: Per-frame rocking-mode rows (``rot`` = φ calc, ``phi_obs``
///     = frame mid-point).
/// :param sigma_m_now_deg: σ_M the latest slices were predicted with.
/// :param osc_deg: Oscillation width of one frame.
/// :param settings: Grouping options (as for :func:`sum_partial_slices`).
/// :param min_reflections: Fewer usable reflections → no estimate.
/// :return: The median σ_M, or ``std::nullopt`` if too few reflections qualify.
[[nodiscard]] std::optional<RockingWidthFit> fit_rocking_width(
    const std::vector<crystal::MosflmReflection>& slices,
    double sigma_m_now_deg,
    double osc_deg,
    const PartialSummationSettings& settings = {},
    std::size_t min_reflections = 40
);

}  // namespace mosflm::diffraction
