#include "mosflm/geometry/MetricTensorProjector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace mosflm {
namespace {

constexpr double kDeg = 3.14159265358979323846 / 180.0;
constexpr double kNullTol = 1e-9;

using Vec6 = Eigen::Matrix<double, 6, 1>;

/// Index pairs for the 6-vector form (g11, g22, g33, g12, g13, g23).
constexpr std::array<std::array<int, 2>, 6> kPairs{{
    {0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2}
}};

Vec6 to_vec6(const Eigen::Matrix3d& G) {
    Vec6 v;
    for (int k = 0; k < 6; ++k) {
        v(k) = G(kPairs[k][0], kPairs[k][1]);
    }
    return v;
}

Eigen::Matrix3d from_vec6(const Vec6& v) {
    Eigen::Matrix3d G;
    for (int k = 0; k < 6; ++k) {
        G(kPairs[k][0], kPairs[k][1]) = v(k);
        G(kPairs[k][1], kPairs[k][0]) = v(k);
    }
    return G;
}

/// Frobenius weights in 6-vector form: off-diagonals appear twice in G.
Eigen::Matrix<double, 6, 6> frobenius_weight() {
    Vec6 w;
    w << 1.0, 1.0, 1.0, 2.0, 2.0, 2.0;
    return w.asDiagonal();
}

/// Kernel of ``G ↦ Rᵀ G R − G`` over all rotations, as orthonormal columns.
Eigen::MatrixXd invariant_basis(const std::vector<Eigen::Matrix3i>& rotations) {
    Eigen::MatrixXd A(6 * std::max<std::size_t>(rotations.size(), 1), 6);
    A.setZero();
    for (std::size_t r = 0; r < rotations.size(); ++r) {
        const Eigen::Matrix3d R = rotations[r].cast<double>();
        for (int k = 0; k < 6; ++k) {
            Vec6 e = Vec6::Zero();
            e(k) = 1.0;
            const Eigen::Matrix3d E = from_vec6(e);
            A.block<6, 1>(static_cast<Eigen::Index>(6 * r), k) =
                to_vec6(R.transpose() * E * R - E);
        }
    }
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
    const auto& s = svd.singularValues();
    const double scale = std::max(1.0, s.size() > 0 ? s(0) : 0.0);
    int rank = 0;
    for (Eigen::Index i = 0; i < s.size(); ++i) {
        rank += s(i) > kNullTol * scale ? 1 : 0;
    }
    return svd.matrixV().rightCols(6 - rank);
}

/// Jacobian of the 6-vector metric with respect to (a,b,c,α,β,γ) in Å and °.
Eigen::Matrix<double, 6, 6> metric_jacobian(const UnitCell& c) {
    const double ca = std::cos(c.alpha_deg * kDeg);
    const double cb = std::cos(c.beta_deg * kDeg);
    const double cg = std::cos(c.gamma_deg * kDeg);
    const double sa = std::sin(c.alpha_deg * kDeg) * kDeg;
    const double sb = std::sin(c.beta_deg * kDeg) * kDeg;
    const double sg = std::sin(c.gamma_deg * kDeg) * kDeg;
    Eigen::Matrix<double, 6, 6> J = Eigen::Matrix<double, 6, 6>::Zero();
    J(0, 0) = 2.0 * c.a;
    J(1, 1) = 2.0 * c.b;
    J(2, 2) = 2.0 * c.c;
    // g12 = a b cosγ
    J(3, 0) = c.b * cg;
    J(3, 1) = c.a * cg;
    J(3, 5) = -c.a * c.b * sg;
    // g13 = a c cosβ
    J(4, 0) = c.c * cb;
    J(4, 2) = c.a * cb;
    J(4, 4) = -c.a * c.c * sb;
    // g23 = b c cosα
    J(5, 1) = c.c * ca;
    J(5, 2) = c.b * ca;
    J(5, 3) = -c.b * c.c * sa;
    return J;
}

std::vector<Eigen::Matrix3i> rotations_of(const gemmi::SpaceGroup* sg) {
    if (sg == nullptr) {
        throw std::invalid_argument("MetricTensorProjector: null space group");
    }
    const gemmi::GroupOps ops = sg->operations();
    std::vector<Eigen::Matrix3i> out;
    out.reserve(ops.sym_ops.size());
    for (const gemmi::Op& op : ops.sym_ops) {
        Eigen::Matrix3i R;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                R(i, j) = op.rot[i][j] / gemmi::Op::DEN;
            }
        }
        out.push_back(R);
    }
    return out;
}

}  // namespace

Eigen::Matrix3d metric_tensor(const UnitCell& cell) {
    const double ca = std::cos(cell.alpha_deg * kDeg);
    const double cb = std::cos(cell.beta_deg * kDeg);
    const double cg = std::cos(cell.gamma_deg * kDeg);
    Eigen::Matrix3d G;
    G << cell.a * cell.a, cell.a * cell.b * cg, cell.a * cell.c * cb,
        cell.a * cell.b * cg, cell.b * cell.b, cell.b * cell.c * ca,
        cell.a * cell.c * cb, cell.b * cell.c * ca, cell.c * cell.c;
    return G;
}

UnitCell cell_from_metric_tensor(const Eigen::Matrix3d& G) {
    const Eigen::LLT<Eigen::Matrix3d> llt(G);
    if (llt.info() != Eigen::Success) {
        throw std::domain_error("metric tensor is not positive definite");
    }
    UnitCell out;
    out.a = std::sqrt(G(0, 0));
    out.b = std::sqrt(G(1, 1));
    out.c = std::sqrt(G(2, 2));
    const auto angle = [](double g, double x, double y) {
        return std::acos(std::clamp(g / (x * y), -1.0, 1.0)) / kDeg;
    };
    out.alpha_deg = angle(G(1, 2), out.b, out.c);
    out.beta_deg = angle(G(0, 2), out.a, out.c);
    out.gamma_deg = angle(G(0, 1), out.a, out.b);
    return out;
}

MetricTensorProjector::MetricTensorProjector(const gemmi::SpaceGroup* sg)
    : MetricTensorProjector(rotations_of(sg)) {}

MetricTensorProjector::MetricTensorProjector(
    const std::vector<Eigen::Matrix3i>& rotations
)
    : basis_(invariant_basis(rotations)) {}

int MetricTensorProjector::free_parameters() const {
    return static_cast<int>(basis_.cols());
}

Eigen::Matrix3d MetricTensorProjector::project(const Eigen::Matrix3d& G) const {
    const Eigen::Matrix<double, 6, 6> W = frobenius_weight();
    const Eigen::MatrixXd BtW = basis_.transpose() * W;
    const Vec6 g = to_vec6(G);
    const Vec6 p = basis_ * (BtW * basis_).ldlt().solve(BtW * g);
    return from_vec6(p);
}

Eigen::VectorXd MetricTensorProjector::coefficients(
    const Eigen::Matrix3d& G
) const {
    const Eigen::Matrix<double, 6, 6> W = frobenius_weight();
    const Eigen::MatrixXd BtW = basis_.transpose() * W;
    return (BtW * basis_).ldlt().solve(BtW * to_vec6(G));
}

Eigen::Matrix3d MetricTensorProjector::metric_from_coefficients(
    const Eigen::VectorXd& c
) const {
    if (c.size() != basis_.cols()) {
        throw std::invalid_argument(
            "metric_from_coefficients: expected " +
            std::to_string(basis_.cols()) + " coefficients"
        );
    }
    return from_vec6(basis_ * c);
}

MetricProjection MetricTensorProjector::project(
    const UnitCell& observed,
    const std::optional<CellSigmas>& sigmas
) const {
    MetricProjection out;
    out.observed = observed;
    out.free_parameters = free_parameters();
    out.dof = 6 - out.free_parameters;

    const Eigen::Matrix3d G = metric_tensor(observed);
    const Vec6 g = to_vec6(G);

    Eigen::Matrix<double, 6, 6> W = frobenius_weight();
    if (sigmas) {
        const std::array<double, 6> s{
            sigmas->a, sigmas->b, sigmas->c,
            sigmas->alpha_deg, sigmas->beta_deg, sigmas->gamma_deg
        };
        if (std::any_of(s.begin(), s.end(), [](double v) { return !(v > 0.0); })) {
            throw std::invalid_argument("cell esds must all be positive");
        }
        Vec6 var;
        for (int k = 0; k < 6; ++k) {
            var(k) = s[k] * s[k];
        }
        const Eigen::Matrix<double, 6, 6> J = metric_jacobian(observed);
        const Eigen::Matrix<double, 6, 6> cov = J * var.asDiagonal() * J.transpose();
        W = cov.ldlt().solve(Eigen::Matrix<double, 6, 6>::Identity());
    }

    const Eigen::MatrixXd BtW = basis_.transpose() * W;
    const Vec6 p = basis_ * (BtW * basis_).ldlt().solve(BtW * g);
    const Eigen::Matrix3d Gc = from_vec6(p);
    out.constrained = cell_from_metric_tensor(Gc);

    const UnitCell& o = observed;
    const UnitCell& c = out.constrained;
    out.max_edge_dev_A = std::max(
        {std::fabs(o.a - c.a), std::fabs(o.b - c.b), std::fabs(o.c - c.c)}
    );
    out.max_angle_dev_deg = std::max(
        {std::fabs(o.alpha_deg - c.alpha_deg),
         std::fabs(o.beta_deg - c.beta_deg),
         std::fabs(o.gamma_deg - c.gamma_deg)}
    );
    out.rel_metric_dev = (G - Gc).norm() / G.norm();

    if (sigmas) {
        const Vec6 r = g - p;
        out.chi2 = r.dot(W * r);
        out.p_value = out.dof > 0 ? chi2_upper_tail(*out.chi2, out.dof) : 1.0;
    }
    return out;
}

double chi2_upper_tail(double chi2, int dof) {
    if (dof < 1) {
        throw std::invalid_argument("chi2_upper_tail: dof must be >= 1");
    }
    if (!(chi2 > 0.0)) {
        return 1.0;
    }
    const double h = 0.5 * chi2;
    // Q(k/2, x/2): even k sums a Poisson tail; odd k starts from erfc.
    double term = 0.0;
    double sum = 0.0;
    if (dof % 2 == 0) {
        term = 1.0;
        for (int i = 0; i < dof / 2; ++i) {
            sum += term;
            term *= h / (i + 1);
        }
        return std::exp(-h) * sum;
    }
    term = std::sqrt(h) / std::tgamma(1.5);
    for (int i = 0; i < (dof - 1) / 2; ++i) {
        sum += term;
        term *= h / (i + 1.5);
    }
    return std::erfc(std::sqrt(h)) + std::exp(-h) * sum;
}

}  // namespace mosflm
