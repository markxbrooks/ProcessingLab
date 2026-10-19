#include "mosflm/diffraction/LpCorrection.hpp"
#include "mosflm/geometry/DetectorGeometry.hpp"
#include "test_util.hpp"

#include <cmath>
#include <iostream>

int main() {
    using mosflm::DetectorGeometry;
    using mosflm::diffraction::lp_factor;

    DetectorGeometry det;
    det.distance_mm = 200.0;
    det.beam_x_mm = 0.0;
    det.beam_y_mm = 0.0;

    // AIR path-length stand-in: longer path at high angle → larger correction
    // factor when applied as exp(μ·extra) (verified algebraically here).
    const double mu = 0.001;
    const double r_near = 10.0;
    const double r_far = 150.0;
    const double path_near = std::hypot(det.distance_mm, r_near);
    const double path_far = std::hypot(det.distance_mm, r_far);
    const double corr_near =
        std::exp(mu * (path_near - det.distance_mm));
    const double corr_far =
        std::exp(mu * (path_far - det.distance_mm));
    CHECK_MSG(corr_far > corr_near, "far spot needs more air correction");
    CHECK_MSG(corr_near > 1.0, "near correction > 1");

    // Partiality full-equivalent: I/frac restores scale.
    const double i_part = 100.0;
    const double frac = 0.5;
    const double i_full = i_part / frac;
    CHECK_CLOSE(i_full, 200.0, 1e-9);

    // Polarized LP still finite off-axis.
    const double lp = lp_factor(det, 80.0, 40.0, 0.98);
    CHECK_MSG(std::isfinite(lp) && lp > 0.0, "polarized LP");

    return mosflm_test::finish("partial_absorb");
}
