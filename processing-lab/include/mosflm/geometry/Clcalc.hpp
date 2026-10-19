#pragma once

#include "mosflm/geometry/UnitCell.hpp"

#include <Eigen/Dense>

namespace mosflm {

/// MOSFLM ``CLCALC``: recover cell parameters from a metric tensor ``G``.
///
/// ``G(i,i) = |axis_i|²``, ``G(i,j) = axis_i · axis_j``. Angles returned in
/// degrees (MOSFLM ``CELL(6)`` order).
[[nodiscard]] UnitCell cell_from_metric_tensor(const Eigen::Matrix3d& G);

/// Build the direct-space metric tensor from a unit cell (inverse of CLCALC).
[[nodiscard]] Eigen::Matrix3d metric_tensor(const UnitCell& cell);

}  // namespace mosflm
