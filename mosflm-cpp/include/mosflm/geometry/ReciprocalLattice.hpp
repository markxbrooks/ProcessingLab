#pragma once

#include "mosflm/geometry/UnitCell.hpp"

#include <Eigen/Dense>

namespace mosflm {

/// Reciprocal-space companion to :class:`UnitCell`.
struct ReciprocalLattice {
    UnitCell cell;  ///< Reciprocal cell (a*,b*,c* in 1/Å; angles in degrees)

    [[nodiscard]] static ReciprocalLattice from_direct(const UnitCell& direct);

    /// Reciprocal basis rows (a*, b*, c*) in the same lab frame as
    /// ``UnitCell::basis_rows()``.
    [[nodiscard]] Eigen::Matrix3d basis_rows() const;

    /// d-spacing for Miller indices (Å).
    [[nodiscard]] double d_spacing(int h, int k, int l) const;
};

}  // namespace mosflm
