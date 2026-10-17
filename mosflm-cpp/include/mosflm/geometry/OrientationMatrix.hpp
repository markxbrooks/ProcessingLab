#pragma once

#include "mosflm/geometry/UnitCell.hpp"

#include <Eigen/Dense>

namespace mosflm {

/// Orientation / UB matrix: maps Miller indices to laboratory reciprocal
/// coordinates (Å⁻¹). Convention: ``r* = UB * h``.
class OrientationMatrix {
public:
    OrientationMatrix() = default;
    explicit OrientationMatrix(Eigen::Matrix3d ub);

    /// Build UB from direct cell + orthogonal crystal rotation R
    /// (``UB = R * A*`` where A* is the reciprocal orthogonalization).
    static OrientationMatrix from_cell_and_rotation(
        const UnitCell& cell,
        const Eigen::Matrix3d& R = Eigen::Matrix3d::Identity()
    );

    [[nodiscard]] const Eigen::Matrix3d& ub() const { return ub_; }
    [[nodiscard]] Eigen::Matrix3d ub_inverse() const;

    /// Laboratory reciprocal vector for (h,k,l).
    [[nodiscard]] Eigen::Vector3d miller_to_lab(int h, int k, int l) const;

    /// Fractional Miller indices for a lab reciprocal vector.
    [[nodiscard]] Eigen::Vector3d lab_to_miller(const Eigen::Vector3d& r) const;

private:
    Eigen::Matrix3d ub_{Eigen::Matrix3d::Identity()};
};

}  // namespace mosflm
