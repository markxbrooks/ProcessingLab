#include "mosflm/diffraction/LpCorrection.hpp"

#include <algorithm>
#include <cmath>

namespace mosflm::diffraction {
namespace {

/// Near-beam |sin 2θ| below this → LP factor returns 1.0 (not physically
/// valid on-axis; avoids singularity).
constexpr double kSinTwoThetaFloor = 1e-8;

}  // namespace

double lp_factor(
    const DetectorGeometry& detector,
    double x_mm,
    double y_mm,
    double polarization_fraction
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
    const double s2sq = s2 * s2;
    const double c2sq = c2 * c2;
    // ρ ∈ [0,1]; 0 → unpolarized.
    const double rho = std::clamp(polarization_fraction, 0.0, 1.0);
    double P = 0.5 * (1.0 + c2sq);
    if (rho > 0.0 && r > 1e-9) {
        // ψ = azimuth in detector plane (horizontal polarization plane).
        const double psi = std::atan2(dy, dx);
        P = 0.5 * ((1.0 + c2sq) - rho * std::cos(2.0 * psi) * s2sq);
        if (P < 1e-6) {
            P = 1e-6;
        }
    }
    const double abs_s2 = std::fabs(s2);
    const double L = (abs_s2 > kSinTwoThetaFloor) ? (1.0 / abs_s2) : 1.0;
    return L * P;
}

void apply_lp_inplace(
    const DetectorGeometry& detector,
    double x_mm,
    double y_mm,
    double& intensity,
    double& sigma,
    double polarization_fraction
) {
    const double lp =
        lp_factor(detector, x_mm, y_mm, polarization_fraction);
    intensity *= lp;
    sigma *= lp;
}

double rotation_lp_correction(
    const Eigen::Vector3d& s0,
    const Eigen::Vector3d& s1,
    const Eigen::Vector3d& axis,
    const Eigen::Vector3d& polarization_normal,
    double polarization_fraction,
    double min_inverse_lorentz
) {
    const double n0 = s0.norm();
    const double n1 = s1.norm();
    const double na = axis.norm();
    const double nn = polarization_normal.norm();
    if (!(n0 > 0.0) || !(n1 > 0.0) || !(na > 0.0) || !(nn > 0.0)) {
        return 1.0;
    }
    const Eigen::Vector3d u0 = s0 / n0;
    const Eigen::Vector3d u1 = s1 / n1;
    const double inv_lorentz = std::max(
        std::fabs((axis / na).dot(u1.cross(u0))), min_inverse_lorentz
    );
    const double p = std::clamp(polarization_fraction, 0.5, 1.0);
    const double cos2t = u0.dot(u1);
    const double ndot = (polarization_normal / nn).dot(u1);
    const double P = std::max(
        (1.0 - 2.0 * p) * (1.0 - ndot * ndot) + p * (1.0 + cos2t * cos2t), 1e-6
    );
    return inv_lorentz / P;
}

}  // namespace mosflm::diffraction
