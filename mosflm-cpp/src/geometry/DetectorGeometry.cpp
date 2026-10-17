#include "mosflm/geometry/DetectorGeometry.hpp"

namespace mosflm {

bool DetectorGeometry::lab_to_detector_mm(
    const Eigen::Vector3d& s1,
    double& x_mm,
    double& y_mm
) const {
    if (std::fabs(s1.z()) < 1e-15) {
        return false;
    }
    const double t = distance_mm / s1.z();
    if (t <= 0.0) {
        return false;
    }
    x_mm = s1.x() * t + beam_x_mm;
    y_mm = s1.y() * t + beam_y_mm;
    return true;
}

void DetectorGeometry::mm_to_pixel(
    double x_mm,
    double y_mm,
    double& xf,
    double& yf
) const {
    xf = x_mm / pixel_size_x_mm;
    yf = y_mm / pixel_size_y_mm;
}

double DetectorGeometry::two_theta_deg(double x_mm, double y_mm) const {
    const double dx = x_mm - beam_x_mm;
    const double dy = y_mm - beam_y_mm;
    const double r = std::sqrt(dx * dx + dy * dy);
    return std::atan2(r, distance_mm) * (180.0 / 3.14159265358979323846);
}

}  // namespace mosflm
