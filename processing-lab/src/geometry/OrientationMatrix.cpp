#include "mosflm/geometry/OrientationMatrix.hpp"

#include "mosflm/geometry/LabFrame.hpp"


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

OrientationMatrix OrientationMatrix::from_mosflm_orientation_matrix(
    const Eigen::Matrix3d& A_mosflm,
    double wavelength_A
) {
    return OrientationMatrix(mosflm_amat_to_cpp_ub(A_mosflm, wavelength_A));
}

UnitCell OrientationMatrix::cell() const {
    // UB = R·A⁻ᵀ, so the rows of UB⁻¹ are the direct axes (R orthogonal).
    // METRIC uses MOSFLM's degree constant, as orthogonalization_matrix.
    return UnitCell::from_basis_rows(ub_.inverse());
}

OrientationMatrix OrientationMatrix::with_cell(const UnitCell& cell) const {
    const Eigen::Matrix3d b_now = from_cell_and_rotation(this->cell()).ub();
    const Eigen::Matrix3d u = ub_ * b_now.inverse();
    // Nearest orthogonal matrix (polar decomposition); W·Vᵀ keeps det(u)'s sign.
    const Eigen::JacobiSVD<Eigen::Matrix3d> svd(u, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Eigen::Matrix3d rot = svd.matrixU() * svd.matrixV().transpose();
    return from_cell_and_rotation(cell, rot);
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
