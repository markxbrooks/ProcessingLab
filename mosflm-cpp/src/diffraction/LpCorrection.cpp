#include "mosflm/diffraction/LpCorrection.hpp"

#include <cmath>

namespace mosflm::diffraction {
namespace {

/// Near-beam |sin 2θ| below this → LP factor returns 1.0 (not physically
/// valid on-axis; avoids singularity for the unpolarized stills stand-in).
constexpr double kSinTwoThetaFloor = 1e-8;

}  // namespace

double lp_factor(
    const DetectorGeometry& detector,
    double x_mm,
    double y_mm
) {
    if (!(detector.distance_mm > 0.0) || !std::isfinite(detector.distance_mm) ||
        !std::isfinite(x_mm) || !std::isfinite(y_mm)) {
        return 1.0;
    }
    const double dx = x_mm - detector.beam_x_mm;
    const double dy = y_mm - detector.beam_y_mm;
    const double r = std::sqrt(dx * dx + dy * dy);
    const double two_theta = std::atan2(r, detector.distance_mm);
    const double c2 = std::cos(two_theta);
    const double s2 = std::sin(two_theta);
    const double P = 0.5 * (1.0 + c2 * c2);  // unpolarized
    const double abs_s2 = std::fabs(s2);
    const double L = (abs_s2 > kSinTwoThetaFloor) ? (1.0 / abs_s2) : 1.0;
    return L * P;
}

void apply_lp_inplace(
    const DetectorGeometry& detector,
    double x_mm,
    double y_mm,
    double& intensity,
    double& sigma
) {
    const double lp = lp_factor(detector, x_mm, y_mm);
    intensity *= lp;
    sigma *= lp;
}

}  // namespace mosflm::diffraction
