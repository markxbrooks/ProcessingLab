#include "mosflm/geometry/UnitCell.hpp"
#include "test_util.hpp"

#include <Eigen/Dense>

using mosflm::UnitCell;

int main() {
    const UnitCell tet{58.0, 58.0, 150.0, 90.0, 90.0, 90.0};

    CHECK_CLOSE(tet.volume(), 58.0 * 58.0 * 150.0, 1e-6);
    CHECK(tet.is_valid());

    const UnitCell recip = tet.reciprocal();
    CHECK_CLOSE(recip.a, 1.0 / 58.0, 1e-12);
    CHECK_CLOSE(recip.b, 1.0 / 58.0, 1e-12);
    CHECK_CLOSE(recip.c, 1.0 / 150.0, 1e-12);
    CHECK_CLOSE(recip.alpha_deg, 90.0, 1e-9);
    CHECK_CLOSE(recip.gamma_deg, 90.0, 1e-9);

    const Eigen::Matrix3d A = tet.orthogonalization_matrix();
    CHECK_CLOSE(A(0, 0), 58.0, 1e-9);
    CHECK_CLOSE(A(1, 1), 58.0, 1e-9);
    CHECK_CLOSE(A(2, 2), 150.0, 1e-9);

    const Eigen::Matrix3d rows = tet.basis_rows();
    CHECK_CLOSE(rows.row(0).norm(), 58.0, 1e-9);
    CHECK_CLOSE(rows.row(2).norm(), 150.0, 1e-9);

    const UnitCell recovered = UnitCell::from_basis_rows(rows);
    CHECK_CLOSE(recovered.a, tet.a, 1e-9);
    CHECK_CLOSE(recovered.b, tet.b, 1e-9);
    CHECK_CLOSE(recovered.c, tet.c, 1e-9);
    CHECK_CLOSE(recovered.alpha_deg, 90.0, 1e-7);
    CHECK_CLOSE(recovered.beta_deg, 90.0, 1e-7);
    CHECK_CLOSE(recovered.gamma_deg, 90.0, 1e-7);
    CHECK_CLOSE(recovered.volume(), tet.volume(), 1e-3);

    // Monoclinic roundtrip exercises the non-orthogonal branch.
    const UnitCell mono{30.0, 40.0, 50.0, 90.0, 105.0, 90.0};
    const UnitCell mono_rt = UnitCell::from_basis_rows(mono.basis_rows());
    CHECK_CLOSE(mono_rt.a, mono.a, 1e-9);
    CHECK_CLOSE(mono_rt.b, mono.b, 1e-9);
    CHECK_CLOSE(mono_rt.c, mono.c, 1e-9);
    CHECK_CLOSE(mono_rt.beta_deg, mono.beta_deg, 1e-7);
    CHECK_CLOSE(mono_rt.volume(), mono.volume(), 1e-3);

    const UnitCell zero_a{0.0, 1.0, 1.0, 90.0, 90.0, 90.0};
    CHECK(!zero_a.is_valid());

    return mosflm_test::finish("unit_cell");
}
