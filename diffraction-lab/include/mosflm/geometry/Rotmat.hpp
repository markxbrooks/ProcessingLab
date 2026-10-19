#pragma once

#include <Eigen/Dense>

namespace mosflm {

/// MOSFLM ``ROTMAT``: ``R = Rz(phiz) · Ry(phiy) · Rx(phix)``.
///
/// Angles in degrees (``degrees == true``) or radians.
[[nodiscard]] Eigen::Matrix3d rotmat_from_missetting(
    double phix,
    double phiy,
    double phiz,
    bool degrees = true
);

}  // namespace mosflm
