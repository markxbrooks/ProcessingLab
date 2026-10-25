#pragma once

#include "mosflm/indexing/SpotList.hpp"
#include "mosflm/indexing/SpotToRlp.hpp"

#include <cstddef>
#include <vector>

namespace mosflm::indexing {

/// Knobs of the DIALS ``find_max_cell`` nearest-neighbour analysis.
struct MaxCellSettings {
    /// Spots are compared only within φ blocks of this width (degrees).
    double step_deg{45.0};
    /// Margin on the most probable neighbour distance (DIALS
    /// ``max_cell_multiplier``).
    double tolerance{1.3};
    /// Histogram bins above this fraction of the tallest bin count as
    /// populated; the last such bin sets the estimate.
    double max_height_fraction{0.25};
    /// Target neighbours per histogram bin.
    int nn_per_bin{5};
    /// Bin ``log10`` of the distances instead (DIALS ``histogram_binning=log``).
    bool log_binning{false};
    /// In (0, 1]: use this percentile of the distances instead of the
    /// histogram rule (DIALS ``nearest_neighbor_percentile``); 0 disables.
    double percentile{0.0};
    /// Drop spots at hexagonal-ice ring resolutions (DIALS ``filter_ice``).
    bool filter_ice{true};
    /// Full width of an ice ring in ``1/d²`` (Å⁻²).
    double ice_width{0.004};
    /// Drop reflections that overlap another on the detector (DIALS
    /// ``filter_overlaps``): a merged chain holding two spots of one frame,
    /// or two reflections within ``overlap_px`` whose φ spans overlap.
    bool filter_overlaps{true};
    double overlap_px{6.0};
    /// Spots this close on the detector (pixels) and in φ (degrees) are
    /// chained into one reflection before the analysis: a MULTIFRAME list
    /// repeats a reflection on each frame it spans.
    double merge_px{3.0};
    double merge_dphi_deg{0.5};
};

/// Result of :func:`estimate_max_cell`.
struct MaxCellEstimate {
    bool ok{false};
    /// Largest expected direct-cell edge (Å), including ``tolerance``.
    double max_cell_A{0.0};
    /// Reflections after merging repeats on adjacent frames.
    std::size_t n_reflections{0};
    /// Reflections rejected at ice-ring resolutions.
    std::size_t n_ice{0};
    /// Reflections rejected as overlapping.
    std::size_t n_overlaps{0};
    /// Neighbour distances that entered the histogram.
    std::size_t n_neighbours{0};
};

/**
 * Hexagonal-ice (P6₃/mmc, a 4.498 Å, c 7.338 Å) reflection positions.
 *
 * \param d_min_A Highest resolution to generate (Å).
 * \return Sorted, distinct ``1/d²`` (Å⁻²) of the allowed reflections.
 */
[[nodiscard]] std::vector<double> ice_ring_d_star_sq(double d_min_A);

/**
 * Largest expected cell edge from reciprocal-lattice nearest neighbours.
 *
 * Port of DIALS ``NeighborAnalysis``: each RLP's nearest neighbour within
 * its group gives a direct-space distance ``1/|Δr|``; the top 1% are
 * dropped, the rest histogrammed, and the far edge of the last bin above
 * ``max_height_fraction`` of the peak, times ``tolerance``, is the estimate.
 * With ``settings.percentile`` set, that percentile of all distances times
 * ``tolerance`` is used instead.
 *
 * \param x,y,z Reciprocal-lattice points (Å⁻¹).
 * \param group Neighbour group of each point (φ block and Ewald side); points
 *        are only compared within a group.
 * \param settings Analysis knobs.
 * \return The estimate; ``ok`` is false with fewer than 10 neighbour pairs.
 */
[[nodiscard]] MaxCellEstimate estimate_max_cell(
    const std::vector<float>& x,
    const std::vector<float>& y,
    const std::vector<float>& z,
    const std::vector<int>& group,
    const MaxCellSettings& settings = {}
);

/**
 * Largest expected cell edge from a ``.spt`` spot list.
 *
 * Port of DIALS ``find_max_cell``. Chains repeats of one reflection on
 * adjacent frames into an intensity-weighted centroid (DIALS works on 3D
 * spots) and maps each to its RLP as AUTOINDEX does. Ice-ring and
 * overlapping reflections are dropped, the rest grouped by
 * ``settings.step_deg`` blocks of φ and by side of the rotation axis
 * (entering / exiting the Ewald sphere), then the nearest-neighbour
 * analysis runs.
 *
 * \param spots Spot list (all spots are used).
 * \param geometry Spot → RLP geometry.
 * \param settings Analysis knobs.
 * \return The estimate.
 */
[[nodiscard]] MaxCellEstimate estimate_max_cell(
    const SpotList& spots,
    const SpotToRlpGeometry& geometry,
    const MaxCellSettings& settings = {}
);

}  // namespace mosflm::indexing
