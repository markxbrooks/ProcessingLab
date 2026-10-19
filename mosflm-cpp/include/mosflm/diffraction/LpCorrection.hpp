#pragma once

#include "mosflm/geometry/DetectorGeometry.hpp"

namespace mosflm::diffraction {

/// Lorentz–polarization helpers for intensity scaling.
///
/// Single source of truth for the C ABI ``mosflm_lp_factor`` and integrators.
///
/// Unpolarized (``polarization_fraction == 0``):
/// ``P = (1 + cos² 2θ) / 2``, ``L = 1 / |sin 2θ|``.
///
/// Synchrotron (Kabsch / XDS ``FRACTION_OF_POLARIZATION`` ρ):
/// ``P = ½ [(1 + cos² 2θ) − ρ cos(2ψ) sin² 2θ]`` where ψ is the azimuth of
/// the diffracted beam in the detector plane relative to the horizontal
/// polarization plane (XDS default).
///
/// When ``|sin 2θ|`` is below an internal floor, returns ``1.0``.
[[nodiscard]] double lp_factor(
    const DetectorGeometry& detector,
    double x_mm,
    double y_mm,
    double polarization_fraction = 0.0
);

/// Scale ``intensity`` / ``sigma`` in place by ``lp_factor``.
void apply_lp_inplace(
    const DetectorGeometry& detector,
    double x_mm,
    double y_mm,
    double& intensity,
    double& sigma,
    double polarization_fraction = 0.0
);

}  // namespace mosflm::diffraction
