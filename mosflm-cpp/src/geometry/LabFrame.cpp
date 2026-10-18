#include "mosflm/geometry/LabFrame.hpp"

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

}  // namespace mosflm
