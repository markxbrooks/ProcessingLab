#pragma once

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/crystal/MtzAdapter.hpp"

#include <iosfwd>
#include <vector>

namespace mosflm::diffraction {

/// One resolution shell of XDS-style completeness / quality statistics.
struct ResolutionShellStats {
    double d_limit_A{0.0};  ///< High-resolution edge of the shell (Å)
    int observed{0};
    int unique{0};
    int possible{0};
    double completeness_pct{0.0};
    double r_obs_pct{-99.9};   ///< Missing → −99.9 (XDS convention)
    double r_exp_pct{-99.9};
    int compared{0};
    double mean_isigma{-99.0};
    double rmeas_pct{-99.9};
    double cc_half_pct{0.0};  ///< CC(1/2) as percentage (−100…100)
    bool cc_significant{false};
    double anomal_corr_pct{0.0};
    double sigano{0.0};
    int nano{0};
};

/// Full completeness / quality report for an unmerged PROCESS reflection list.
struct DataQualityReport {
    std::vector<ResolutionShellStats> shells;
    ResolutionShellStats total{};
    int n_in_subset{0};          ///< Observations with I/σ ≥ min_isigma
    int n_rejected_misfits{0};   ///< flag ≠ 0 (MOSFLM misfit bits)
    int n_systematic_absent{0};
    /// Mean I/σ of systematically absent observations (CORRECT-lite sanity).
    double mean_isigma_systematic_absent{-99.0};
    int n_accepted_obs{0};
    int n_unique_accepted{0};
    /// The rest of the input rows, by the filter that dropped them, so that
    /// rows = rejected_misfits + invalid_sigma + systematic_absent +
    /// below_min_isigma + outside_resolution + accepted_obs.
    int n_invalid_sigma{0};       ///< σ ≤ 0, or I or σ not finite
    int n_below_min_isigma{0};    ///< I/σ < ``min_isigma``
    int n_outside_resolution{0};  ///< d outside [``dmin_A``, ``dmax_A``]
    /// The range and cut actually used (defaults replace an empty range).
    double dmax_A{0.0};
    double dmin_A{0.0};
    double min_isigma{0.0};
};

/// R-meas contributions by image, using the same filters and denominator as
/// the cumulative PROCESS data-quality report.
struct ImageRmeasStats {
    int image{0};
    int n_accepted{0};
    int n_compared{0};
    double rmeas_pct{-99.9};
};

/// Build XDS-style completeness / R-factor / R-meas / CC(1/2) table from
/// unmerged MOSFLM PROCESS rows.
///
/// Anomalous columns are left at zero (FRIEDEL'S_LAW assumed; no anomalous
/// pairing in PROCESS). Shells are equal-width in ``1/d`` between
/// ``dmax_A`` and ``dmin_A``.
[[nodiscard]] DataQualityReport compute_data_quality(
    const crystal::CrystalMetadata& crystal,
    const std::vector<crystal::MosflmReflection>& reflections,
    double dmax_A,
    double dmin_A,
    int n_shells = 9,
    double min_isigma = -3.0
);

/// Compute per-image R-meas against each reflection's mean over the complete
/// input list. This makes image values additive into the cumulative R-meas.
[[nodiscard]] std::vector<ImageRmeasStats> compute_image_rmeas(
    const crystal::CrystalMetadata& crystal,
    const std::vector<crystal::MosflmReflection>& reflections,
    double dmax_A,
    double dmin_A,
    double min_isigma = -3.0
);

/// Print the XDS ``COMPLETENESS AND QUALITY OF DATA SET`` block (legend +
/// resolution table + observation counts).
///
/// The numeric table sits between ``RESOLUTION_BINS`` and
/// ``RESOLUTION_BINS_END``. The marker line names the columns. ``#`` lines
/// are the aligned header and notes. Each data line is one shell, then a
/// ``total`` line. Fields are space-separated and carry no ``%`` or ``*``;
/// ``cc_sig`` is 1 when CC(1/2) is significant at the 0.1% level, and
/// R-factors are −99.9 when ``compared`` is 0.
void write_data_quality_table(
    std::ostream& os,
    const DataQualityReport& report
);

}  // namespace mosflm::diffraction
