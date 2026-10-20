#pragma once

#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/image/DiffractionImage.hpp"

#include <cstdint>
#include <vector>

namespace mosflm::diffraction {

/// INIT/DEFPIX-lite: circular + resolution trusted-pixel map for fine-pixel
/// detectors (Eiger). Untrusted pixels are excluded from background rims.
struct TrustedRegionSettings {
    /// Keep pixels inside this fraction of the panel half-diagonal (0–1).
    /// Ignored when ``cover_panel`` is true.
    double trusted_radius_frac{0.98};
    /// Trust the rectangular panel. The beam-centered circle misses the far
    /// corner when the beam is offset. ``dmin_A`` / ``dmax_A`` still apply.
    bool cover_panel{false};
    /// Resolution limits in Å (dmax, dmin); 0 disables that bound.
    double dmax_A{0.0};
    double dmin_A{0.0};
    /// Beam centre in mm (panel origin); if ≤0 use detector.beam_* or panel mid.
    double beam_x_mm{0.0};
    double beam_y_mm{0.0};
};

/// ``trusted[row * width + col]`` — 1 = use, 0 = reject.
[[nodiscard]] std::vector<std::uint8_t> build_trusted_mask(
    const image::DiffractionImage& image,
    const DetectorGeometry& detector,
    const TrustedRegionSettings& settings
);

[[nodiscard]] inline bool is_trusted(
    const std::vector<std::uint8_t>& mask,
    int width,
    int col,
    int row
) {
    if (mask.empty() || width <= 0) {
        return true;
    }
    const int h = static_cast<int>(mask.size()) / width;
    if (col < 0 || row < 0 || col >= width || row >= h) {
        return false;
    }
    return mask[static_cast<std::size_t>(row) * static_cast<std::size_t>(width) +
                static_cast<std::size_t>(col)] != 0;
}

}  // namespace mosflm::diffraction
