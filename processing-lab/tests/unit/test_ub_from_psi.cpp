#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/geometry/UnitCell.hpp"
#include "mosflm/legacy/MosflmCApi.h"
#include "test_util.hpp"

using mosflm::OrientationMatrix;
using mosflm::UnitCell;

int main() {
    mosflm_cell6 cell{};
    cell.a = 58.0;
    cell.b = 58.0;
    cell.c = 150.0;
    cell.alpha_deg = 90.0;
    cell.beta_deg = 90.0;
    cell.gamma_deg = 90.0;

    // Zero PSI == canonical R = I orientation.
    const double zero[3] = {0.0, 0.0, 0.0};
    mosflm_ub9 ub0{};
    CHECK_EQ(mosflm_ub_from_psi(&cell, zero, &ub0), MOSFLM_SUCCESS);

    const UnitCell uc{58.0, 58.0, 150.0, 90.0, 90.0, 90.0};
    const Eigen::Matrix3d canon =
        OrientationMatrix::from_cell_and_rotation(uc).ub();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            CHECK_CLOSE(ub0.m[3 * i + j], canon(i, j), 1e-9);
        }
    }

    // Non-zero PSI changes the orientation.
    const double psi[3] = {0.0, 0.0, 90.0};
    mosflm_ub9 ubz{};
    CHECK_EQ(mosflm_ub_from_psi(&cell, psi, &ubz), MOSFLM_SUCCESS);
    double delta = 0.0;
    for (int i = 0; i < 9; ++i) {
        delta += std::fabs(ubz.m[i] - ub0.m[i]);
    }
    CHECK(delta > 0.05);

    // Explicit matrix copy (row-major).
    const double mat[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    mosflm_ub9 ubm{};
    CHECK_EQ(mosflm_ub_from_matrix9(mat, &ubm), MOSFLM_SUCCESS);
    for (int i = 0; i < 9; ++i) {
        CHECK_CLOSE(ubm.m[i], mat[i], 1e-12);
    }

    // NULL arguments rejected.
    CHECK_EQ(mosflm_ub_from_psi(nullptr, zero, &ub0),
             MOSFLM_ERROR_INVALID_ARGUMENT);
    CHECK_EQ(mosflm_ub_from_matrix9(nullptr, &ubm),
             MOSFLM_ERROR_INVALID_ARGUMENT);

    return mosflm_test::finish("ub_from_psi");
}
