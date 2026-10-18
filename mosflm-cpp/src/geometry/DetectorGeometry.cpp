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
    // Beam-relative flat projection in the C++ lab frame (beam along +z).
    double xr = s1.x() * t;
    double yr = s1.y() * t;

    // First-order TILT/TWIST in the MOSFLM XG/YG frame (XG↔y_mm, YG↔x_mm).
    // Full ROTMAT2 detector normal is later Phase D work.
    const double tilt = tilt_mrad * 1.0e-3;
    const double twist = twist_mrad * 1.0e-3;
    double xg = yr;
    double yg = xr;
    const double xg_out = xg + twist * yg;
    const double yg_out = yg + tilt * xg;

    x_mm = yg_out + beam_x_mm;
    y_mm = xg_out + beam_y_mm;
    return true;
}

void DetectorGeometry::mm_to_pixel(
    double x_mm,
    double y_mm,
    double& xf,
    double& yf
) const {
    // MOSFLM YSCALE stretches YG (↔ C++ x_mm): YPIX = Y_mm * YSCAL / RAST
    // so the effective pitch on that axis is pixel_size / yscale.
    const double px =
        (std::fabs(yscale) > 1e-15) ? (pixel_size_x_mm / yscale)
                                    : pixel_size_x_mm;
    xf = x_mm / px;
    yf = y_mm / pixel_size_y_mm;
}

double DetectorGeometry::two_theta_deg(double x_mm, double y_mm) const {
    const double dx = x_mm - beam_x_mm;
    const double dy = y_mm - beam_y_mm;
    const double r = std::sqrt(dx * dx + dy * dy);
    return std::atan2(r, distance_mm) * (180.0 / 3.14159265358979323846);
}

}  // namespace mosflm
