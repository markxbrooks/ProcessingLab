#include "mosflm/geometry/Dstar.hpp"
#include "mosflm/geometry/Rotmat.hpp"
#include "test_util.hpp"

#include <Eigen/Dense>

using mosflm::Dstar;
using mosflm::OrientationMatrix;
using mosflm::UnitCell;
using mosflm::rotmat_from_missetting;

int main() {
    const UnitCell tet{58.0, 58.0, 150.0, 90.0, 90.0, 90.0};
    const OrientationMatrix om = OrientationMatrix::from_cell_and_rotation(tet);

    const Eigen::Vector3d q = Dstar::q_vector(om, 1, 0, 0);
    CHECK_CLOSE(q.x(), 1.0 / 58.0, 1e-12);
    CHECK_CLOSE(Dstar::d_spacing(q), 58.0, 1e-9);

    CHECK_CLOSE(Dstar::d_spacing(tet, 1, 0, 0), 58.0, 1e-9);
    CHECK_CLOSE(Dstar::d_spacing(tet, 0, 1, 0), 58.0, 1e-9);
    CHECK_CLOSE(Dstar::d_spacing(tet, 0, 0, 1), 150.0, 1e-9);

    // On the Ewald sphere the excitation error is ~0 for a reflection whose
    // |q| equals the sphere chord; just check the sign/order of magnitude.
    const Eigen::Vector3d qz = Dstar::q_vector(om, 0, 0, 1);
    const double err = Dstar::excitation_error(qz, 1.0);
    CHECK(std::isfinite(err));

    // Rotation matrices.
    const Eigen::Matrix3d I = rotmat_from_missetting(0.0, 0.0, 0.0);
    CHECK_CLOSE(I(0, 0), 1.0, 1e-12);
    CHECK_CLOSE(I(1, 1), 1.0, 1e-12);
    CHECK_CLOSE(I(2, 2), 1.0, 1e-12);

    // Rz(90 deg): x -> y.
    const Eigen::Matrix3d rz = rotmat_from_missetting(0.0, 0.0, 90.0);
    const Eigen::Vector3d vx = rz * Eigen::Vector3d(1.0, 0.0, 0.0);
    CHECK_CLOSE(vx.x(), 0.0, 1e-8);
    CHECK_CLOSE(vx.y(), 1.0, 1e-12);

    // Rx(90 deg): y -> z.
    const Eigen::Matrix3d rx = rotmat_from_missetting(90.0, 0.0, 0.0);
    const Eigen::Vector3d vy = rx * Eigen::Vector3d(0.0, 1.0, 0.0);
    CHECK_CLOSE(vy.y(), 0.0, 1e-8);
    CHECK_CLOSE(vy.z(), 1.0, 1e-12);

    // Radian input matches degree input.
    const Eigen::Matrix3d rz_rad =
        rotmat_from_missetting(0.0, 0.0, 3.14159265358979323846 / 2.0, false);
    CHECK_CLOSE(rz_rad(0, 0), rz(0, 0), 1e-8);
    CHECK_CLOSE(rz_rad(0, 1), rz(0, 1), 1e-8);

    return mosflm_test::finish("dstar_rotmat");
}
