#pragma once

#include "mosflm/geometry/UnitCell.hpp"

#include <Eigen/Dense>
#include <gemmi/symmetry.hpp>

#include <optional>
#include <vector>

namespace mosflm {

/// Standard uncertainties of the six cell parameters (Å and degrees).
struct CellSigmas {
    double a{};
    double b{};
    double c{};
    double alpha_deg{};
    double beta_deg{};
    double gamma_deg{};
};

/// Direct-space metric tensor ``G_ij = a_i · a_j`` (Å²) of a unit cell.
[[nodiscard]] Eigen::Matrix3d metric_tensor(const UnitCell& cell);

/// Recover cell parameters from a metric tensor.
///
/// \throws std::domain_error if ``G`` is not symmetric positive definite.
[[nodiscard]] UnitCell cell_from_metric_tensor(const Eigen::Matrix3d& G);

/// Observed cell versus the best cell allowed by a candidate's symmetry.
///
/// Kept separate from the DPS / SUGGLATT penalty: it measures only how far
/// the observed metric departs from the candidate lattice metric.
struct MetricProjection {
    UnitCell observed;
    UnitCell constrained;
    int free_parameters{6};    ///< Dimension of the invariant metric subspace.
    double max_edge_dev_A{};   ///< max |Δa|, |Δb|, |Δc| (Å).
    double max_angle_dev_deg{};///< max |Δα|, |Δβ|, |Δγ| (°).
    double rel_metric_dev{};   ///< ‖G_obs − G_con‖_F / ‖G_obs‖_F.
    /// Weighted χ² of the metric misfit; set only when esds were supplied.
    std::optional<double> chi2;
    int dof{0};                ///< ``6 − free_parameters``.
    /// Upper-tail probability P(χ²_dof ≥ chi2); set with ``chi2``.
    std::optional<double> p_value;
};

/// Projects a measured cell onto the metric allowed by a symmetry group.
///
/// The allowed metrics are the symmetric ``G`` with ``Rᵀ G R = G`` for every
/// rotational part ``R`` (acting on fractional coordinates) of the group.
/// That set is a linear subspace of the six independent ``G`` components; it
/// is found numerically from the operations, so no per-crystal-system special
/// cases are needed. Monoclinic unique axis and trigonal/rhombohedral setting
/// follow from the operations of the space group as supplied.
///
/// The constrained metric is the least-squares projection of the observed
/// metric onto that subspace, weighted by the inverse covariance propagated
/// from cell esds when given, else by the Frobenius inner product.
class MetricTensorProjector {
public:
    /// Build from the rotations of a Gemmi space group.
    ///
    /// \throws std::invalid_argument if ``sg`` is null.
    explicit MetricTensorProjector(const gemmi::SpaceGroup* sg);

    /// Build from integer rotation matrices on fractional coordinates.
    explicit MetricTensorProjector(const std::vector<Eigen::Matrix3i>& rotations);

    /// Dimension of the invariant metric subspace (1 cubic … 6 triclinic).
    [[nodiscard]] int free_parameters() const;

    /// Project a metric tensor onto the invariant subspace (Frobenius).
    [[nodiscard]] Eigen::Matrix3d project(const Eigen::Matrix3d& G) const;

    /// Coordinates of the Frobenius projection of ``G`` in the invariant
    /// basis (length :meth:`free_parameters`). Use with
    /// :meth:`metric_from_coefficients` to refine only symmetry-allowed
    /// metric directions.
    [[nodiscard]] Eigen::VectorXd coefficients(const Eigen::Matrix3d& G) const;

    /// Metric tensor for basis coordinates ``c``.
    ///
    /// \throws std::invalid_argument if ``c`` has the wrong length.
    [[nodiscard]] Eigen::Matrix3d metric_from_coefficients(
        const Eigen::VectorXd& c
    ) const;

    /// Constrain a cell and report its departure from the allowed metric.
    ///
    /// \param observed Measured cell in the group's conventional setting.
    /// \param sigmas   Optional cell esds; all must be positive.
    /// \throws std::invalid_argument if any supplied esd is not positive.
    /// \throws std::domain_error if the projected metric is not positive
    ///     definite (observed cell far from this lattice type).
    [[nodiscard]] MetricProjection project(
        const UnitCell& observed,
        const std::optional<CellSigmas>& sigmas = std::nullopt
    ) const;

private:
    /// Columns span the invariant subspace in 6-vector form
    /// ``(g11, g22, g33, g12, g13, g23)``.
    Eigen::MatrixXd basis_;
};

/// Upper-tail χ² probability for integer ``dof`` ≥ 1 (closed form).
[[nodiscard]] double chi2_upper_tail(double chi2, int dof);

}  // namespace mosflm
