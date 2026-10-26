#pragma once

#include "mosflm/image/DiffractionImage.hpp"
#include "mosflm/indexing/SpotList.hpp"

#include <cstdint>
#include <memory>
#include <vector>

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
    /// Pixels the search may use (non-zero), row-major, e.g. DEFPIX's
    /// :meth:`DefpixResult::spot_search_mask`. Masked pixels count neither
    /// towards the background statistics nor as peaks. Ignored unless it
    /// matches the image size.
    std::shared_ptr<const std::vector<std::uint8_t>> mask;
    /// Local background (counts per frame), row-major; ``< 0`` where the
    /// global threshold applies. Where set, e.g. DEFPIX's beam halo, a peak
    /// must also clear ``background + threshold_sigma·√(gain·background)``,
    /// and the pixel stays out of the global background statistics: the
    /// flare's noise is not taken for spots, yet real low-resolution spots
    /// on it are kept. Ignored unless it matches the image size.
    std::shared_ptr<const std::vector<float>> local_background;
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
