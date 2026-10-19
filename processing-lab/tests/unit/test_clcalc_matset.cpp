#include "mosflm/geometry/Clcalc.hpp"
#include "mosflm/geometry/Matset.hpp"
#include "test_util.hpp"

#include <Eigen/Dense>

using mosflm::cell_from_metric_tensor;
using mosflm::metric_tensor;
using mosflm::matset;
using mosflm::orientation_from_matset;
using mosflm::UnitCell;

int main() {
    const UnitCell tet{58.0, 58.0, 150.0, 90.0, 90.0, 90.0};
    const Eigen::Matrix3d G = metric_tensor(tet);
    CHECK_CLOSE(G(0, 0), 58.0 * 58.0, 1e-6);
    CHECK_CLOSE(G(2, 2), 150.0 * 150.0, 1e-6);
    CHECK_CLOSE(G(0, 1), 0.0, 1e-4);

    const UnitCell back = cell_from_metric_tensor(G);
    CHECK_CLOSE(back.a, tet.a, 1e-7);
    CHECK_CLOSE(back.b, tet.b, 1e-7);
    CHECK_CLOSE(back.c, tet.c, 1e-7);
    CHECK_CLOSE(back.alpha_deg, 90.0, 1e-6);
    CHECK_CLOSE(back.beta_deg, 90.0, 1e-6);
    CHECK_CLOSE(back.gamma_deg, 90.0, 1e-6);

    // Triclinic roundtrip exercises the general angle branch.
    const UnitCell tri{30.0, 40.0, 50.0, 80.0, 95.0, 110.0};
    const UnitCell tri_rt = cell_from_metric_tensor(metric_tensor(tri));
    CHECK_CLOSE(tri_rt.a, tri.a, 1e-7);
    CHECK_CLOSE(tri_rt.b, tri.b, 1e-7);
    CHECK_CLOSE(tri_rt.c, tri.c, 1e-7);
    CHECK_CLOSE(tri_rt.alpha_deg, tri.alpha_deg, 1e-6);
    CHECK_CLOSE(tri_rt.beta_deg, tri.beta_deg, 1e-6);
    CHECK_CLOSE(tri_rt.gamma_deg, tri.gamma_deg, 1e-6);

    // MATSET: identity U gives UB = B = reciprocal orthogonalization.
    const auto ms = matset(tet);
    CHECK_CLOSE(ms.A(0, 0), 1.0 / 58.0, 1e-12);
    CHECK_CLOSE(ms.A(2, 2), 1.0 / 150.0, 1e-12);
    CHECK_CLOSE(ms.U(0, 0), 1.0, 1e-12);

    // orientation_from_matset matches from_cell_and_rotation.
    const auto om = orientation_from_matset(tet);
    const auto om_direct = mosflm::OrientationMatrix::from_cell_and_rotation(tet);
    CHECK_CLOSE(om.ub()(0, 0), om_direct.ub()(0, 0), 1e-12);
    CHECK_CLOSE(om.ub()(2, 2), om_direct.ub()(2, 2), 1e-12);

    return mosflm_test::finish("clcalc_matset");
}
