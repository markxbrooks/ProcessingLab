#include "mosflm/geometry/Clcalc.hpp"

#include <cmath>
#include <stdexcept>

namespace mosflm {
namespace {

// MOSFLM trigparam.f
constexpr double kDtor = 0.0174532925;
constexpr double kRtod = 57.2957795131;

double clamp_cos(double x) {
    if (x > 1.0) {
        return 1.0;
    }
    if (x < -1.0) {
        return -1.0;
    }
    return x;
}

}  // namespace

UnitCell cell_from_metric_tensor(const Eigen::Matrix3d& G) {
    // Mirrors SUBROUTINE CLCALC(CELL, A) with A = metric tensor.
    UnitCell cell;
    for (int i = 0; i < 3; ++i) {
        if (G(i, i) <= 0.0) {
            throw std::runtime_error(
                "cell_from_metric_tensor: non-positive G(i,i)"
            );
        }
    }
    cell.a = std::sqrt(G(0, 0));
    cell.b = std::sqrt(G(1, 1));
    cell.c = std::sqrt(G(2, 2));
    cell.alpha_deg = kRtod * std::acos(clamp_cos(G(1, 2) / (cell.b * cell.c)));
    cell.beta_deg = kRtod * std::acos(clamp_cos(G(0, 2) / (cell.a * cell.c)));
    cell.gamma_deg = kRtod * std::acos(clamp_cos(G(0, 1) / (cell.a * cell.b)));
    (void)kDtor;
    return cell;
}

Eigen::Matrix3d metric_tensor(const UnitCell& cell) {
    const double ca = std::cos(cell.alpha_deg * kDtor);
    const double cb = std::cos(cell.beta_deg * kDtor);
    const double cg = std::cos(cell.gamma_deg * kDtor);
    Eigen::Matrix3d G = Eigen::Matrix3d::Zero();
    G(0, 0) = cell.a * cell.a;
    G(1, 1) = cell.b * cell.b;
    G(2, 2) = cell.c * cell.c;
    G(0, 1) = G(1, 0) = cell.a * cell.b * cg;
    G(0, 2) = G(2, 0) = cell.a * cell.c * cb;
    G(1, 2) = G(2, 1) = cell.b * cell.c * ca;
    return G;
}

}  // namespace mosflm
