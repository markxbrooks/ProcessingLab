#pragma once

#include "mosflm/geometry/DetectorGeometry.hpp"

namespace mosflm::diffraction {

/// Lorentz–polarization helpers for intensity scaling (unpolarized beam).
///
/// Single source of truth for the C ABI ``mosflm_lp_factor`` and integrators.
///
/// ``P = (1 + cos² 2θ) / 2``, ``L = 1 / |sin 2θ|``. When ``|sin 2θ|`` is below
/// an internal floor, returns ``1.0`` (avoids singularity; **not** a physically
/// valid on-axis correction). Combined ``LP = L · P`` multiplies intensities
/// when ``apply_lp`` is enabled on the integrate / profile path.
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
