#pragma once

#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/geometry/UnitCell.hpp"

#include <Eigen/Dense>

namespace mosflm {

/// Compact reciprocal-space helpers used on the PREDICT path (Fortran ``DSTAR``
/// family: |d*| and related quantities without the full COMMON-heavy routine).
struct Dstar {
    /// Reciprocal-space vector ``q = UB · (h,k,l)`` (Å⁻¹).
    [[nodiscard]] static Eigen::Vector3d q_vector(
        const OrientationMatrix& ub,
        int h,
        int k,
        int l
    );

    /// d-spacing (Å) from ``|q|``.
    [[nodiscard]] static double d_spacing(const Eigen::Vector3d& q);

    /// d-spacing from cell + Miller indices via reciprocal orthogonalization.
    [[nodiscard]] static double d_spacing(
        const UnitCell& cell,
        int h,
        int k,
        int l
    );

    /// Ewald excitation error ``|s0 + q| - 1/λ`` (Å⁻¹) with ``s0 = (0,0,1/λ)``.
    [[nodiscard]] static double excitation_error(
        const Eigen::Vector3d& q,
        double wavelength_A
    );
};

}  // namespace mosflm
