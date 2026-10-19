#pragma once

#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/geometry/UnitCell.hpp"

#include <Eigen/Dense>

namespace mosflm {

/// Result of a thin ``MATSET``-style construction (no COMMON / no derivatives).
struct MatsetResult {
    Eigen::Matrix3d B{Eigen::Matrix3d::Identity()};  ///< Reciprocal orth. (Å⁻¹)
    Eigen::Matrix3d U{Eigen::Matrix3d::Identity()};  ///< Crystal orientation
    Eigen::Matrix3d A{Eigen::Matrix3d::Identity()};  ///< ``A = U · B`` (lab UB)
};

/// Build ``B`` (reciprocal orthogonalization) and ``A = U · B``.
///
/// Full Fortran ``MATSET`` also fills ``DUBDCC`` cell derivatives into COMMON;
/// those are deferred. This covers the PREDICT-path need: UB from cell + U.
[[nodiscard]] MatsetResult matset(
    const UnitCell& cell,
    const Eigen::Matrix3d& U = Eigen::Matrix3d::Identity()
);

/// Convenience: wrap :func:`matset` as an :class:`OrientationMatrix`.
[[nodiscard]] OrientationMatrix orientation_from_matset(
    const UnitCell& cell,
    const Eigen::Matrix3d& U = Eigen::Matrix3d::Identity()
);

}  // namespace mosflm
