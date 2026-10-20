#pragma once

#include "mosflm/indexing/SpotList.hpp"

#include <cstddef>
#include <vector>

namespace mosflm::indexing {

/// Geometry needed to map ``.spt`` spots → reciprocal lab vectors (Å⁻¹).
struct SpotToRlpGeometry {
    double distance_mm{160.0};
    double wavelength_A{0.9790};
    double beam_x_mm{94.0};
    double beam_y_mm{94.7};
};

/**
 * MOSFLM DPSCHOOSE-style reciprocal coordinates from spot centroids.
 *
 * FR_MAT (omega + invflag) → triangular Ewald map → ROTMATY about Y using
 * ``phi_end_deg`` (matches AutoIndexer / INDEXFIT xyz).
 *
 * \param keep Indices into ``spots.spots`` (empty → use all spots).
 * \param[out] x,y,z Reciprocal lab components (Å⁻¹), same length.
 * \return Number of reciprocal vectors written.
 */
[[nodiscard]] std::size_t spots_to_rlp_xyz(
    const SpotList& spots,
    const std::vector<std::size_t>& keep,
    const SpotToRlpGeometry& geometry,
    std::vector<float>& x,
    std::vector<float>& y,
    std::vector<float>& z
);

}  // namespace mosflm::indexing
