#pragma once

#include "mosflm/image/DiffractionImage.hpp"
#include "mosflm/indexing/SpotList.hpp"

namespace mosflm::indexing {

/// Settings for the lightweight DiffractionLab-style peak finder.
struct SpotFindSettings {
    /// Peak threshold in robust background σ (MAD × 1.4826).
    float threshold_sigma{5.0f};
    /// Minimum separation between accepted peaks (pixels).
    int min_distance{8};
    /// Cap on returned spots (strongest first). Default sized for
    /// AUTOINDEX on large cells (single-frame ≤500 often yields sublattices).
    int max_spots{1500};
    /// Border (pixels) excluded from search.
    int border{20};
    /// φ start / end (degrees). If ``phi_end <= phi_start`` and the image has
    /// an oscillation range, that range is used.
    double phi_start_deg{0.0};
    double phi_end_deg{0.0};
    /// INIT-lite: subtract estimated blank (ADU) before peak search.
    double blank_adu{0.0};
    /// INIT-lite: scale background σ by ``sqrt(gain)`` (XDS GAIN.cbf stand-in).
    double gain{1.0};
};

/// DiffractionLab-equivalent spot finder (MAD background + 3×3 local maxima).
///
/// No OpenCV dependency. Suitable for SMV / CBF / Eiger frames already loaded
/// as :class:`image::DiffractionImage`. Writes MOSFLM ``.spt`` via
/// :meth:`SpotList::write_spt`.
class SpotFinder {
public:
    /// Find spots; coordinates are millimetres from the image origin.
    [[nodiscard]] SpotList find(
        const image::DiffractionImage& image,
        const SpotFindSettings& settings = {}
    ) const;
};

}  // namespace mosflm::indexing
