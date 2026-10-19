#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <string>

namespace mosflm {

/// Direct-space unit cell (Å and degrees), matching MOSFLM CELL(6) ordering.
struct UnitCell {
    double a{1.0};
    double b{1.0};
    double c{1.0};
    double alpha_deg{90.0};
    double beta_deg{90.0};
    double gamma_deg{90.0};

    /// Cell volume in Å³ (MOSFLM ``INVCEL`` formula). Returns ≤0 on failure.
    [[nodiscard]] double volume() const;

    /// Reciprocal cell parameters (1/Å and degrees), MOSFLM ``INVCEL``.
    [[nodiscard]] UnitCell reciprocal() const;

    /// Orthogonalization matrix A such that fractional → Cartesian (Å):
    /// ``x = A * f``. Conventional crystallographic setting with a along x,
    /// b in the xy plane.
    [[nodiscard]] Eigen::Matrix3d orthogonalization_matrix() const;

    /// Basis vectors stored **row-wise** (MOSFLM ``METRIC`` convention):
    /// row 0 = a⃗, row 1 = b⃗, row 2 = c⃗ in an orthonormal lab frame.
    [[nodiscard]] Eigen::Matrix3d basis_rows() const;

    /// Recover cell parameters from row-wise basis (MOSFLM ``METRIC``).
    [[nodiscard]] static UnitCell from_basis_rows(const Eigen::Matrix3d& B);

    [[nodiscard]] bool is_valid(double eps = 1e-8) const;
    [[nodiscard]] std::string describe() const;
};

}  // namespace mosflm
