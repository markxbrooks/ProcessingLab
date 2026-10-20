#include "mosflm/geometry/DetectorGeometry.hpp"

#include "mosflm/geometry/LabFrame.hpp"

#include <algorithm>

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

double detector_corner_dmin_A(
    const DetectorGeometry& det,
    double beam_x_mm,
    double beam_y_mm,
    double margin_px
) {
    const double px = det.pixel_size_x_mm;
    const double py = det.pixel_size_y_mm;
    const double dist = det.distance_mm;
    const double wave = det.wavelength_A;
    if (!(px > 0.0) || !(py > 0.0) || !(dist > 0.0) || !(wave > 0.0) ||
        det.nfast < 2 || det.nslow < 2) {
        return 0.0;
    }
    const double m = std::max(0.0, margin_px);
    const double c0 = m;
    const double c1 = static_cast<double>(det.nfast - 1) - m;
    const double r0 = m;
    const double r1 = static_cast<double>(det.nslow - 1) - m;
    if (!(c1 > c0) || !(r1 > r0)) {
        return 0.0;
    }

    auto d_at = [&](double col, double row) -> double {
        const double rad = std::hypot(col * px - beam_x_mm, row * py - beam_y_mm);
        if (rad < 1e-9) {
            return 1.0e6;
        }
        const double theta = 0.5 * std::atan2(rad, dist);
        const double s = std::sin(theta);
        if (s < 1e-12) {
            return 1.0e6;
        }
        return wave / (2.0 * s);
    };

    double dmin = d_at(c0, r0);
    dmin = std::min(dmin, d_at(c1, r0));
    dmin = std::min(dmin, d_at(c0, r1));
    dmin = std::min(dmin, d_at(c1, r1));
    // Keep a reflection whose d equals the corner pixel (`d < dmin` rejects).
    return dmin * (1.0 - 1.0e-4);
}

bool spot_on_panel(
    const DetectorGeometry& det,
    double x_mm,
    double y_mm,
    double margin_px,
    bool use_panel_beam,
    double panel_beam_x_mm,
    double panel_beam_y_mm
) {
    const double margin = std::max(0.0, margin_px);
    if (use_panel_beam) {
        double xdet = 0.0;
        double ydet = 0.0;
        cpp_mm_to_mosflm_xdet_yscale(
            x_mm,
            y_mm,
            panel_beam_x_mm,
            panel_beam_y_mm,
            det.pixel_size_x_mm,
            det.pixel_size_y_mm,
            det.yscale,
            xdet,
            ydet
        );
        const double fast_hi = static_cast<double>(det.nfast) - margin;
        const double slow_hi = static_cast<double>(det.nslow) - margin;
        return xdet >= margin && xdet < fast_hi && ydet >= margin &&
               ydet < slow_hi;
    }
    double xf = 0.0;
    double yf = 0.0;
    det.mm_to_pixel(x_mm, y_mm, xf, yf);
    const double half_f = 0.5 * static_cast<double>(det.nfast) - margin;
    const double half_s = 0.5 * static_cast<double>(det.nslow) - margin;
    return std::fabs(xf) < half_f && std::fabs(yf) < half_s;
}

}  // namespace mosflm
