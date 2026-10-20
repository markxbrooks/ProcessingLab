/**
 * Index-parity supercell reduction: 3× axis (thaumatin trap), 2× face
 * diagonal, and a no-op on the true cell.
 */

#include "mosflm/indexing/SupercellReduce.hpp"
#include "test_util.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <random>
#include <vector>

using mosflm::indexing::reduce_supercell_by_index_parity;

namespace {

void fill_rlps(
    const Eigen::Matrix3d& B,
    std::vector<float>& x,
    std::vector<float>& y,
    std::vector<float>& z
) {
    std::mt19937 rng(11);
    std::uniform_int_distribution<int> hkl(-8, 8);
    std::normal_distribution<double> noise(0.0, 0.00003);
    for (int n = 0; n < 600; ++n) {
        const Eigen::Vector3d h(hkl(rng), hkl(rng), hkl(rng));
        if (h.squaredNorm() < 0.5) {
            continue;
        }
        const Eigen::Vector3d r = B * h;
        x.push_back(static_cast<float>(r.x() + noise(rng)));
        y.push_back(static_cast<float>(r.y() + noise(rng)));
        z.push_back(static_cast<float>(r.z() + noise(rng)));
    }
}

}  // namespace

int main() {
    // Thaumatin P41212 57.9 × 57.9 × 151 (8HVE); x = B · hkl.
    Eigen::Matrix3d B = Eigen::Matrix3d::Zero();
    B(0, 0) = 1.0 / 57.9;
    B(1, 1) = 1.0 / 57.9;
    B(2, 2) = 1.0 / 151.0;
    const double det_true = std::fabs(B.determinant());
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> z;
    fill_rlps(B, x, y, z);

    // Supercell AMAT = B · M⁻¹ (M integer, det M = n): indices become M·hkl.
    auto run = [&](const Eigen::Matrix3d& M, int expect_factor,
                   const char* name) {
        const Eigen::Matrix3d A_super = B * M.inverse();
        const auto sr = reduce_supercell_by_index_parity(A_super, x, y, z);
        CHECK_MSG(sr.factor == expect_factor,
                  std::string(name) + " factor=" + mosflm_test::fmt(sr.factor));
        CHECK_MSG(
            std::fabs(std::fabs(sr.amat.determinant()) / det_true - 1.0) < 1e-6,
            std::string(name) + " volume not restored"
        );
    };

    // 58 × 151 × 174 (b tripled; axes permuted as in the DPS table).
    Eigen::Matrix3d M3;
    M3 << 1, 0, 0,
          0, 0, 1,
          0, 3, 0;
    run(M3, 3, "axis3");

    // √2 basal: (a+b), (a−b), c → index 2.
    Eigen::Matrix3d M2;
    M2 << 1, 1, 0,
          1, -1, 0,
          0, 0, 1;
    run(M2, 2, "diag2");

    // 6× (2 along a, 3 along c) sheds both.
    Eigen::Matrix3d M6;
    M6 << 2, 0, 0,
          0, 1, 0,
          0, 0, 3;
    run(M6, 6, "axis6");

    run(Eigen::Matrix3d::Identity(), 1, "true");

    return mosflm_test::finish("supercell_reduce");
}
