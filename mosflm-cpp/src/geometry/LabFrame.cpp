#include "mosflm/geometry/LabFrame.hpp"

namespace mosflm {

Eigen::Matrix3d mosflm_lab_to_cpp() {
    Eigen::Matrix3d R;
    // Columns: images of MOSFLM e_x, e_y, e_z in the C++ lab frame.
    // +X (beam / photon travel) → +Z, +Y → −Y, +Z (spindle) → +X.
    R << 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0;
    return R;
}

Eigen::Matrix3d mosflm_amat_to_cpp_ub(const Eigen::Matrix3d& A_mosflm) {
    return mosflm_lab_to_cpp() * A_mosflm;
}

}  // namespace mosflm
