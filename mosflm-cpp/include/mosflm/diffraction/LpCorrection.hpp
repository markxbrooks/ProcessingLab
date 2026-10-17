#pragma once

#include "mosflm/geometry/DetectorGeometry.hpp"

namespace mosflm::diffraction {

/// Lorentz–polarization helpers for intensity scaling (unpolarized beam).
///
/// ``P = (1 + cos² 2θ) / 2``, ``L = 1 / |sin 2θ|`` (clamped near the beam).
/// Combined factor ``LP = L · P`` multiplies observed intensities when
/// converting to a roughly corrected scale (stand-in for MOSFLM LP).
[[nodiscard]] double lp_factor(
    const DetectorGeometry& detector,
    double x_mm,
    double y_mm
);

/// Scale ``intensity`` / ``sigma`` in place by ``lp_factor``.
void apply_lp_inplace(
    const DetectorGeometry& detector,
    double x_mm,
    double y_mm,
    double& intensity,
    double& sigma
);

}  // namespace mosflm::diffraction
