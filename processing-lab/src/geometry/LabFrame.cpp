#include "mosflm/geometry/LabFrame.hpp"

#include <cmath>

namespace mosflm {

Eigen::Matrix3d mosflm_lab_to_cpp() {
    Eigen::Matrix3d R;
    // Columns: images of MOSFLM e_x, e_y, e_z in the C++ lab frame.
    // +X (beam / photon travel) → +Z, +Y → −Y, +Z (spindle) → +X.
    R << 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0;
    return R;
}

Eigen::Matrix3d mosflm_amat_to_cpp_ub(
    const Eigen::Matrix3d& A_mosflm,
    double wavelength_A
) {
    if (wavelength_A <= 0.0) {
        return mosflm_lab_to_cpp() * A_mosflm;
    }
    return mosflm_lab_to_cpp() * (A_mosflm / wavelength_A);
}

Eigen::Matrix3d cpp_ub_to_mosflm_amat(
    const Eigen::Matrix3d& ub_cpp,
    double wavelength_A
) {
    const Eigen::Matrix3d Rinv = mosflm_lab_to_cpp().inverse();
    if (wavelength_A <= 0.0) {
        return Rinv * ub_cpp;
    }
    return Rinv * ub_cpp * wavelength_A;
}

Eigen::Matrix3d dps_xyz_to_cpp_lab(bool invertx) {
    Eigen::Matrix3d Q;
    if (invertx) {
        Q << 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    } else {
        Q << 0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    }
    return Q;
}

double ccom_from_spt_omega(double omega_deg, bool invertx) {
    double c = invertx ? -omega_deg : (omega_deg + 180.0);
    while (c < 0.0) {
        c += 360.0;
    }
    while (c >= 360.0) {
        c -= 360.0;
    }
    return (std::fabs(c) < 1e-12) ? 0.0 : c;  // no "-0" in logs
}

Eigen::Matrix3d xyz_reciprocal_to_mosflm_amat(
    const Eigen::Matrix3d& om_t_xyz,
    double wavelength_A,
    bool invertx
) {
    return cpp_ub_to_mosflm_amat(
        dps_xyz_to_cpp_lab(invertx) * om_t_xyz, wavelength_A
    );
}

Eigen::Matrix3d mosflm_amat_to_xyz_reciprocal(
    const Eigen::Matrix3d& A_mosflm,
    double wavelength_A,
    bool invertx
) {
    return dps_xyz_to_cpp_lab(invertx).inverse() *
           mosflm_amat_to_cpp_ub(A_mosflm, wavelength_A);
}

void cpp_mm_to_mosflm_xdet(
    double x_mm,
    double y_mm,
    double beam_x_mm,
    double beam_y_mm,
    double pixel_size_x_mm,
    double pixel_size_y_mm,
    double& xdet_px,
    double& ydet_px
) {
    // XG ↔ +y_mm, YG ↔ +x_mm (see LabFrame.hpp).
    xdet_px = beam_x_mm / pixel_size_x_mm + y_mm / pixel_size_x_mm;
    ydet_px = beam_y_mm / pixel_size_y_mm + x_mm / pixel_size_y_mm;
}

void cpp_mm_to_mosflm_xdet_yscale(
    double x_mm,
    double y_mm,
    double beam_x_mm,
    double beam_y_mm,
    double pixel_size_x_mm,
    double pixel_size_y_mm,
    double yscale,
    double& xdet_px,
    double& ydet_px
) {
    const double py =
        (std::fabs(yscale) > 1e-15) ? (pixel_size_y_mm / yscale)
                                    : pixel_size_y_mm;
    cpp_mm_to_mosflm_xdet(
        x_mm, y_mm, beam_x_mm, beam_y_mm, pixel_size_x_mm, py, xdet_px, ydet_px
    );
}

}  // namespace mosflm
