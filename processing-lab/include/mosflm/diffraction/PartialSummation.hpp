#pragma once

#include "mosflm/crystal/MtzAdapter.hpp"

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

}  // namespace mosflm::diffraction
