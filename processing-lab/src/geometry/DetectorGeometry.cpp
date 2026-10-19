#include "mosflm/geometry/DetectorGeometry.hpp"

namespace mosflm {
namespace {

// MOSFLM trigparam.f DTOR (Clemens Vonrhein).
constexpr double kDtor = 0.0174532925;

}  // namespace

Eigen::Vector3d rotmat2_detnor(double tilt_rad, double twist_rad) {
    // ROTMAT2: C1=cos(-TWIST), C2=cos(-TILT), S1=sin(-TWIST), S2=sin(-TILT)
    const double c1 = std::cos(-twist_rad);
    const double c2 = std::cos(-tilt_rad);
    const double s1 = std::sin(-twist_rad);
    const double s2 = std::sin(-tilt_rad);
    return {s2, s1 * c2, c1 * c2};
}

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
    const double xr = s1.x() * t;
    const double yr = s1.y() * t;

    // MOSFLM XG/YG (ideal, at distance_mm): XG↔y_mm, YG↔x_mm.
    const double xg = yr;
    const double yg = xr;

    // MMTOPX works in 10 µm units; XTOFD = distance_mm * 100.
    const double xtofd = distance_mm * 100.0;
    if (xtofd < 1e-6) {
        return false;
    }

    double dxtof = xtofra;
    if (use_rotmat2) {
        // NEWT / NUTWIST: DXTOF = 1 / (XRED·n1 + YRED·n2 + n3)
        const double tilt_rad = tilt_mrad * 1.0e-3;
        const double twist_rad = twist_mrad * 1.0e-3;
        const Eigen::Vector3d n = rotmat2_detnor(tilt_rad, twist_rad);
        const double xred = (xg * 100.0) / xtofd;
        const double yred = (yg * 100.0) / xtofd;
        const double denom = xred * n.x() + yred * n.y() + n.z();
        if (std::fabs(denom) < 1e-15) {
            return false;
        }
        dxtof = 1.0 / denom;
    } else if (std::fabs(tilt_mrad) > 1e-15 || std::fabs(twist_mrad) > 1e-15) {
        // Old MMTOPX (default): TILT = ITILT * FDIST, FDIST = DTOR/(XTOFD*100)
        const double fdist = kDtor / (xtofd * 100.0);
        const double tilt = tilt_mrad * fdist;
        const double twist = twist_mrad * fdist;
        const double X = xg * 100.0;
        const double Y = yg * 100.0;
        dxtof = xtofra + tilt * X + twist * Y;
    }

    // CCOM (OMEGA0) rotation of ideal XG/YG before distance scale.
    const double om = ccom_deg * kDtor;
    const double cos_om = std::cos(om);
    const double sin_om = std::sin(om);
    const double xm = xg * cos_om - yg * sin_om;
    const double ym = xg * sin_om + yg * cos_om;

    const double xg_out = dxtof * xm;
    const double yg_out = dxtof * ym;

    // Back to C++ beam-relative mm (YG↔x_mm, XG↔y_mm), then beam origin.
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
