#include "mosflm/geometry/OrientationMatrix.hpp"
#include "test_util.hpp"

#include <Eigen/Dense>

using mosflm::OrientationMatrix;
using mosflm::UnitCell;

int main() {
    const UnitCell tet{58.0, 58.0, 150.0, 90.0, 90.0, 90.0};
    const OrientationMatrix om =
        OrientationMatrix::from_cell_and_rotation(tet);

    // Reciprocal basis vector a* = 1/a along x.
    const Eigen::Vector3d r100 = om.miller_to_lab(1, 0, 0);
    CHECK_CLOSE(r100.x(), 1.0 / 58.0, 1e-12);
    CHECK_CLOSE(r100.y(), 0.0, 1e-8);
    CHECK_CLOSE(r100.z(), 0.0, 1e-8);

    const Eigen::Vector3d r001 = om.miller_to_lab(0, 0, 1);
    CHECK_CLOSE(r001.z(), 1.0 / 150.0, 1e-12);

    // Roundtrip: Miller -> lab -> Miller.
    const Eigen::Vector3d hkl(3.0, -2.0, 5.0);
    const Eigen::Vector3d lab = om.ub() * hkl;
    const Eigen::Vector3d back = om.lab_to_miller(lab);
    CHECK_CLOSE(back.x(), hkl.x(), 1e-9);
    CHECK_CLOSE(back.y(), hkl.y(), 1e-9);
    CHECK_CLOSE(back.z(), hkl.z(), 1e-9);

    // ub_inverse * ub = I.
    const Eigen::Matrix3d prod = om.ub_inverse() * om.ub();
    CHECK_CLOSE(prod(0, 0), 1.0, 1e-9);
    CHECK_CLOSE(prod(1, 1), 1.0, 1e-9);
    CHECK_CLOSE(prod(2, 2), 1.0, 1e-9);
    CHECK_CLOSE(prod(0, 1), 0.0, 1e-9);
    CHECK_CLOSE(prod(1, 2), 0.0, 1e-9);

    // A 90 degree rotation about z should swap the x/y reciprocal axes.
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
    R(0, 0) = 0.0;
    R(0, 1) = -1.0;
    R(1, 0) = 1.0;
    R(1, 1) = 0.0;
    const OrientationMatrix rot =
        OrientationMatrix::from_cell_and_rotation(tet, R);
    const Eigen::Vector3d ra = rot.miller_to_lab(1, 0, 0);
    CHECK_CLOSE(ra.x(), 0.0, 1e-8);
    CHECK_CLOSE(ra.y(), 1.0 / 58.0, 1e-12);

    // cell(): the cell a UB implies; with_cell(): same orientation, new cell.
    {
        const UnitCell mono{60.0, 70.0, 80.0, 90.0, 105.0, 90.0};
        const Eigen::Matrix3d R =
            (Eigen::AngleAxisd(0.4, Eigen::Vector3d(1, 2, 3).normalized())).toRotationMatrix();
        const auto ub = OrientationMatrix::from_cell_and_rotation(mono, R);
        const UnitCell c = ub.cell();
        CHECK_CLOSE(c.a, 60.0, 1e-9);
        CHECK_CLOSE(c.b, 70.0, 1e-9);
        CHECK_CLOSE(c.c, 80.0, 1e-9);
        CHECK_CLOSE(c.beta_deg, 105.0, 1e-9);
        CHECK_CLOSE(c.alpha_deg, 90.0, 1e-9);

        const UnitCell other{61.0, 69.5, 81.0, 90.0, 104.0, 90.0};
        const auto moved = ub.with_cell(other);
        const auto want = OrientationMatrix::from_cell_and_rotation(other, R);
        CHECK((moved.ub() - want.ub()).norm() < 1e-12);
        CHECK_CLOSE(moved.cell().c, 81.0, 1e-9);

        // Improper (ProcessingLab's XDS-frame handedness) stays improper.
        Eigen::Matrix3d flip = R;
        flip.row(0) *= -1.0;
        const auto improper = OrientationMatrix::from_cell_and_rotation(mono, flip);
        const auto improper_moved = improper.with_cell(other);
        CHECK(improper_moved.ub().determinant() < 0.0);
        CHECK((improper_moved.ub() - OrientationMatrix::from_cell_and_rotation(other, flip).ub()).norm() < 1e-12);
    }

    return mosflm_test::finish("orientation");
}
