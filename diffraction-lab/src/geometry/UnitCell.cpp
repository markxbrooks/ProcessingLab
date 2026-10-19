#include "mosflm/geometry/UnitCell.hpp"

#include <cstdio>

namespace mosflm {
namespace {

constexpr double kDeg2Rad = 1.0 / 57.29578;  // MOSFLM constant
constexpr double kRad2Deg = 57.29578;

double axis_length(const UnitCell& cell, int idx) {
    if (idx == 0) {
        return cell.a;
    }
    if (idx == 1) {
        return cell.b;
    }
    return cell.c;
}

}  // namespace

double UnitCell::volume() const {
    const double ca = std::cos(alpha_deg * kDeg2Rad);
    const double cb = std::cos(beta_deg * kDeg2Rad);
    const double cc = std::cos(gamma_deg * kDeg2Rad);
    const double arg = 1.0 + 2.0 * ca * cb * cc - ca * ca - cb * cb - cc * cc;
    if (arg <= 0.0) {
        return 0.0;
    }
    return a * b * c * std::sqrt(arg);
}

UnitCell UnitCell::reciprocal() const {
    // Mirrors SUBROUTINE INVCEL in mosflm/torefix_noX.f
    UnitCell r{};
    const double ca[3] = {
        std::cos(alpha_deg * kDeg2Rad),
        std::cos(beta_deg * kDeg2Rad),
        std::cos(gamma_deg * kDeg2Rad),
    };
    const double sa[3] = {
        std::sin(alpha_deg * kDeg2Rad),
        std::sin(beta_deg * kDeg2Rad),
        std::sin(gamma_deg * kDeg2Rad),
    };
    const double cell[3] = {a, b, c};
    const double V = volume();
    if (V <= 0.0) {
        throw std::runtime_error("UnitCell::reciprocal: invalid cell volume");
    }
    for (int i = 0; i < 3; ++i) {
        const int j = (i + 1) % 3;
        const int k = (i + 2) % 3;
        const double ri = cell[j] * cell[k] * sa[i] / V;
        double arg =
            kRad2Deg * std::acos((ca[j] * ca[k] - ca[i]) / (sa[j] * sa[k]));
        const int rounded = static_cast<int>(std::lround(arg));
        if (std::fabs(arg - rounded) < 0.0009) {
            arg = static_cast<double>(rounded);
        }
        if (i == 0) {
            r.a = ri;
            r.alpha_deg = arg;
        } else if (i == 1) {
            r.b = ri;
            r.beta_deg = arg;
        } else {
            r.c = ri;
            r.gamma_deg = arg;
        }
    }
    return r;
}

Eigen::Matrix3d UnitCell::orthogonalization_matrix() const {
    const double cos_a = std::cos(alpha_deg * kDeg2Rad);
    const double cos_b = std::cos(beta_deg * kDeg2Rad);
    const double cos_g = std::cos(gamma_deg * kDeg2Rad);
    const double sin_g = std::sin(gamma_deg * kDeg2Rad);
    const double V = volume();
    if (V <= 0.0 || std::fabs(sin_g) < 1e-15) {
        throw std::runtime_error(
            "UnitCell::orthogonalization_matrix: invalid cell"
        );
    }
    Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
    A(0, 0) = a;
    A(0, 1) = b * cos_g;
    A(0, 2) = c * cos_b;
    A(1, 0) = 0.0;
    A(1, 1) = b * sin_g;
    A(1, 2) = c * (cos_a - cos_b * cos_g) / sin_g;
    A(2, 0) = 0.0;
    A(2, 1) = 0.0;
    A(2, 2) = V / (a * b * sin_g);
    return A;
}

Eigen::Matrix3d UnitCell::basis_rows() const {
    return orthogonalization_matrix().transpose();
}

UnitCell UnitCell::from_basis_rows(const Eigen::Matrix3d& B) {
    // Mirrors SUBROUTINE METRIC in mosflm/torefix_noX.f
    UnitCell cell{};
    cell.a = B.row(0).norm();
    cell.b = B.row(1).norm();
    cell.c = B.row(2).norm();

    auto angle_deg = [&](int i, int j) -> double {
        const double r = axis_length(cell, i) * axis_length(cell, j);
        if (r <= 0.0) {
            throw std::runtime_error(
                "UnitCell::from_basis_rows: zero-length axis"
            );
        }
        double cosang = B.row(i).dot(B.row(j)) / r;
        if (std::fabs(cosang) >= 1.0) {
            throw std::runtime_error(
                "UnitCell::from_basis_rows: illegal basis"
            );
        }
        return kRad2Deg * std::acos(cosang);
    };

    // METRIC: (I,J)=(2,1)→gamma; (3,1)→beta; (3,2)→alpha  (1-based)
    cell.gamma_deg = angle_deg(1, 0);
    cell.beta_deg = angle_deg(2, 0);
    cell.alpha_deg = angle_deg(2, 1);
    return cell;
}

bool UnitCell::is_valid(double eps) const {
    return volume() > eps && a > eps && b > eps && c > eps;
}

std::string UnitCell::describe() const {
    char buf[160];
    std::snprintf(
        buf,
        sizeof(buf),
        "a=%.4f b=%.4f c=%.4f α=%.2f β=%.2f γ=%.2f V=%.2f",
        a,
        b,
        c,
        alpha_deg,
        beta_deg,
        gamma_deg,
        volume()
    );
    return std::string(buf);
}

}  // namespace mosflm
