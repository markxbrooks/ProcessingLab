/**
 * Rotation-method Lorentz-polarization correction (DIALS / XDS form).
 */
#include "mosflm/diffraction/LpCorrection.hpp"
#include "test_util.hpp"

#include <cmath>

int main() {
    using mosflm::diffraction::rotation_lp_correction;
    using mosflm_test::close;
    const double kPi = 3.14159265358979323846;
    const Eigen::Vector3d s0(0.0, 0.0, 1.0);
    const Eigen::Vector3d axis(1.0, 0.0, 0.0);
    const Eigen::Vector3d normal(0.0, 1.0, 0.0);

    // Scattering in the vertical (y-z) plane, perpendicular to the axis:
    // 1/L = sin 2θ; normal·ŝ1 = sin 2θ so P = 1 for p = 1.
    const double tt = 20.0 * kPi / 180.0;
    const Eigen::Vector3d s1v(0.0, std::sin(tt), std::cos(tt));
    CHECK(close(rotation_lp_correction(s0, s1v, axis, normal, 1.0), std::sin(tt), 1e-12));
    // Unpolarized: P = (1 + cos² 2θ)/2.
    const double punp = 0.5 * (1.0 + std::cos(tt) * std::cos(tt));
    CHECK(close(rotation_lp_correction(s0, s1v, axis, normal, 0.5), std::sin(tt) / punp, 1e-12));
    CHECK(close(rotation_lp_correction(s0, s1v, axis, normal, 0.0), std::sin(tt) / punp, 1e-12));

    // Scattering in the horizontal (x-z) plane, containing the axis:
    // 1/L floors; polarized P = cos² 2θ.
    const Eigen::Vector3d s1h(std::sin(tt), 0.0, std::cos(tt));
    CHECK(close(
        rotation_lp_correction(s0, s1h, axis, normal, 1.0, 1e-3),
        1e-3 / (std::cos(tt) * std::cos(tt)), 1e-12
    ));

    // Low-angle reflections get a *smaller* correction than high-angle ones.
    const Eigen::Vector3d s1lo(0.0, std::sin(0.02), std::cos(0.02));
    CHECK(rotation_lp_correction(s0, s1lo, axis, normal, 0.98) <
          rotation_lp_correction(s0, s1v, axis, normal, 0.98));
    // Scale-invariant in |s0|, |s1|.
    CHECK(close(
        rotation_lp_correction(2.0 * s0, 3.0 * s1v, axis, normal, 0.98),
        rotation_lp_correction(s0, s1v, axis, normal, 0.98), 1e-12
    ));
    return mosflm_test::finish("rotation_lp");
}
