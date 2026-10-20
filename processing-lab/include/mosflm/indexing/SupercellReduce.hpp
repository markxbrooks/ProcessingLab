#pragma once

#include <Eigen/Dense>

#include <vector>

namespace mosflm::indexing {

/// Index-parity supercell test (DIALS/Mosflm-style sublattice check).
///
/// A basis that is an n× supercell of the true lattice indexes every spot,
/// but the indices obey a hidden congruence m·hkl ≡ 0 (mod n) — e.g. all
/// k ≡ 0 (mod 3) for a cell tripled along b (thaumatin 58×149×174 vs the
/// true 58×58×151). Fit-fraction scoring cannot see this; the parity can.
struct SupercellReduceSettings {
    /// Max |hkl − round(hkl)| for a spot to vote.
    double hkldev{0.30};
    /// Fraction of voting spots that must obey the congruence.
    double min_frac{0.97};
    /// Minimum voting spots before any reduction is attempted.
    int min_indexed{30};
};

struct SupercellReduceResult {
    /// MOSFLM-style AMAT of the reduced basis (x = A · hkl).
    Eigen::Matrix3d amat{Eigen::Matrix3d::Identity()};
    /// Integer basis change: ``amat = amat_in · U`` (det U = ``factor``).
    Eigen::Matrix3d U{Eigen::Matrix3d::Identity()};
    /// Volume reduction factor (1 = unchanged).
    int factor{1};
};

/// Shed supercell multiples of ``amat`` (x = A · hkl) whose index parity
/// shows the true lattice is a sublattice of index 2 or 3. Repeats until no
/// congruence holds. The result is not Niggli-reduced.
[[nodiscard]] SupercellReduceResult reduce_supercell_by_index_parity(
    const Eigen::Matrix3d& amat,
    const std::vector<float>& x,
    const std::vector<float>& y,
    const std::vector<float>& z,
    const SupercellReduceSettings& settings = {}
);

}  // namespace mosflm::indexing
