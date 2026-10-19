#include "mosflm/geometry/ReciprocalLattice.hpp"

namespace mosflm {

ReciprocalLattice ReciprocalLattice::from_direct(const UnitCell& direct) {
    ReciprocalLattice rl;
    rl.cell = direct.reciprocal();
    return rl;
}

Eigen::Matrix3d ReciprocalLattice::basis_rows() const {
    return cell.basis_rows();
}

double ReciprocalLattice::d_spacing(int h, int k, int l) const {
    const Eigen::Vector3d r =
        cell.orthogonalization_matrix() * Eigen::Vector3d(h, k, l);
    const double len = r.norm();
    if (len <= 0.0) {
        return 0.0;
    }
    return 1.0 / len;
}

}  // namespace mosflm
