#include "mosflm/geometry/Matset.hpp"

namespace mosflm {

MatsetResult matset(const UnitCell& cell, const Eigen::Matrix3d& U) {
    // B maps Miller → reciprocal Cartesian in the crystal frame.
    // Matches OrientationMatrix::from_cell_and_rotation with R = I:
    //   A_direct^{-T} = reciprocal orthogonalization.
    const Eigen::Matrix3d A_direct = cell.orthogonalization_matrix();
    const Eigen::Matrix3d B = A_direct.inverse().transpose();
    MatsetResult out;
    out.B = B;
    out.U = U;
    out.A = U * B;
    return out;
}

OrientationMatrix orientation_from_matset(
    const UnitCell& cell,
    const Eigen::Matrix3d& U
) {
    return OrientationMatrix(matset(cell, U).A);
}

}  // namespace mosflm
