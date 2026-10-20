#pragma once

#include "mosflm/geometry/DetectorGeometry.hpp"

#include <Eigen/Dense>

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

/// Rotation-method Lorentz-polarization *correction* (multiply raw I by it).
///
/// Returns ``1 / (L·P)`` with the rotation Lorentz factor
/// ``1/L = |m·(s1 × s0)| / (|s0||s1|)`` and the Kabsch polarization factor
/// ``P = (1 − 2p)(1 − (n·ŝ1)²) + p(1 + cos² 2θ)``, as in DIALS / XDS.
///
/// :param s0: Incident wavevector (any length).
/// :param s1: Diffracted wavevector (any length).
/// :param axis: Rotation axis direction.
/// :param polarization_normal: Normal to the polarization plane.
/// :param polarization_fraction: Fraction ``p`` of the beam polarized in the
///     plane (0.5 = unpolarized; values below 0.5 are treated as 0.5).
/// :param min_inverse_lorentz: Floor on ``1/L`` near the rotation axis.
/// :return: Correction factor; ``I_corrected = I_raw · factor``.
[[nodiscard]] double rotation_lp_correction(
    const Eigen::Vector3d& s0,
    const Eigen::Vector3d& s1,
    const Eigen::Vector3d& axis,
    const Eigen::Vector3d& polarization_normal,
    double polarization_fraction,
    double min_inverse_lorentz = 1e-3
);

}  // namespace mosflm::diffraction
