#pragma once

#include "mosflm/indexing/SpotFinder.hpp"
#include "mosflm/indexing/SpotList.hpp"

#include <filesystem>
#include <functional>
#include <vector>

namespace mosflm::indexing {

/// COLSPOT-lite: accumulate 2D spots over a frame range and merge nearby
/// centroids (3D adjacency approximated by φ proximity + xy distance).
struct MultiFrameSpotSettings {
    SpotFindSettings per_frame{};
    /// Merge spots if xy distance (mm) below this and |Δφ| below phi_tol_deg.
    double merge_xy_mm{0.3};
    double merge_phi_tol_deg{0.5};
    /// Minimum number of frames a merged spot must appear on.
    int min_frames{1};
    int max_spots{2000};
};

/// Load frames via ``read_frame(frame_1based) → DiffractionImage``.
using FrameLoader = std::function<image::DiffractionImage(int)>;

[[nodiscard]] SpotList find_spots_multiframe(
    int first_frame,
    int last_frame,
    const FrameLoader& load,
    const MultiFrameSpotSettings& settings = {}
);

}  // namespace mosflm::indexing
