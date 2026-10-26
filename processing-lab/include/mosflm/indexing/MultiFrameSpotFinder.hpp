#pragma once

#include "mosflm/indexing/SpotFinder.hpp"
#include "mosflm/indexing/SpotList.hpp"

#include <filesystem>
#include <functional>
#include <vector>

namespace mosflm::indexing {

/// COLSPOT-lite: find spots on each frame of a range and merge each
/// reflection's frames into one 3D spot (DIALS-style), so indexing sees one
/// position and one φ per reflection.
struct MultiFrameSpotSettings {
    SpotFindSettings per_frame{};
    /// Track radius (mm): a spot continues a reflection's track from the
    /// previous frame only if it lies within this distance of the track's
    /// last spot (and within four of its running centroid). Must be finite and
    /// ``≥ 0``; ``0`` keeps every per-frame spot, each at its frame's
    /// centre φ. A tolerance, not proof of identity: matches are also
    /// one-to-one and unambiguous (see :func:`find_spots_multiframe`).
    double merge_xy_mm{0.0};
    /// No longer used: tracks follow a reflection across any number of
    /// frames (kept for CONTROL ``FINDSPOTS MERGE_PHI``).
    double merge_phi_tol_deg{0.5};
    /// Minimum number of frames a merged spot must span (``≥ 1``).
    int min_frames{1};
    /// Spots returned: merged reflections when merging, per-frame spots
    /// otherwise, brightest first. ``≤ 0`` means no limit.
    int max_spots{2000};
    /// OpenMP threads for per-frame load+find (``< 1`` is taken as 1).
    int nproc{1};
};

/// Load frames via ``read_frame(frame_1based) → DiffractionImage``. With
/// ``nproc > 1`` it is called from several threads at once, so it must be
/// safe for concurrent calls (the HDF5 reader serialises itself; CBF and
/// SMV frames are separate files).
using FrameLoader = std::function<image::DiffractionImage(int)>;

/**
 * Find spots on frames ``first_frame``…``last_frame`` and, with
 * ``merge_xy_mm > 0``, merge each reflection's frames into one spot.
 *
 * Tracks are built frame by frame. A spot extends a track from the previous
 * frame only when it is that track's sole candidate within ``merge_xy_mm``
 * and the track is the spot's sole candidate (no shared or ambiguous
 * matches), and it also lies within four ``merge_xy_mm`` of the track's
 * running centroid (so a chain of small steps cannot run on indefinitely,
 * while a reflection's real drift over its rocking curve is followed). A track
 * takes at most one spot per frame; a frame without a match ends it, and a
 * later sighting starts a new track.
 *
 * A merged spot is a representative point for indexing, not an integrated
 * reflection: its position and φ are intensity-weighted means over the
 * track (φ from each frame's centre; φ is unwrapped, as in a rotation
 * series, so no 0°/360° folding is done), ``phi_start_deg`` and
 * ``phi_end_deg`` both hold that φ, and its intensity and σ are the
 * brightest frame's, so ``AUTOINDEX THRESHOLD`` keeps its per-frame meaning.
 * Spots are ordered brightest first, ties broken by position, so the output
 * does not depend on thread scheduling. The header (detector, pixel size,
 * beam) is the first frame's; it describes the whole series.
 *
 * :raises std::invalid_argument: If ``last_frame < first_frame``,
 *     ``merge_xy_mm`` is negative or not finite, or ``min_frames < 1``.
 */
[[nodiscard]] SpotList find_spots_multiframe(
    int first_frame,
    int last_frame,
    const FrameLoader& load,
    const MultiFrameSpotSettings& settings = {}
);

}  // namespace mosflm::indexing
