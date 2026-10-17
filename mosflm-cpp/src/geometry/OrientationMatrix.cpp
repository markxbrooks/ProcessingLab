#include "mosflm/geometry/OrientationMatrix.hpp"

namespace mosflm {

OrientationMatrix::OrientationMatrix(Eigen::Matrix3d ub) : ub_(std::move(ub)) {}

OrientationMatrix OrientationMatrix::from_cell_and_rotation(
    const UnitCell& cell,
    const Eigen::Matrix3d& R
) {
    // A maps fractional → Cartesian; A^{-T} maps Miller → reciprocal Cartesian.
    const Eigen::Matrix3d A = cell.orthogonalization_matrix();
    const Eigen::Matrix3d Astar = A.inverse().transpose();
    return OrientationMatrix(R * Astar);
}

Eigen::Matrix3d OrientationMatrix::ub_inverse() const {
    return ub_.inverse();
}

Eigen::Vector3d OrientationMatrix::miller_to_lab(int h, int k, int l) const {
    return ub_ * Eigen::Vector3d(h, k, l);
}

Eigen::Vector3d OrientationMatrix::lab_to_miller(const Eigen::Vector3d& r) const {
    return ub_inverse() * r;
}

}  // namespace mosflm
