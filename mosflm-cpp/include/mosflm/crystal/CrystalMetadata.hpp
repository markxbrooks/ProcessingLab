#pragma once

#include "mosflm/geometry/UnitCell.hpp"

#include <gemmi/symmetry.hpp>
#include <gemmi/unitcell.hpp>

#include <stdexcept>
#include <string>

namespace mosflm::crystal {

/// Crystallographic metadata backed by Gemmi; adapters to mosflm::UnitCell.
struct CrystalMetadata {
    gemmi::UnitCell cell;
    const gemmi::SpaceGroup* spacegroup{nullptr};

    [[nodiscard]] static CrystalMetadata from_parameters(
        double a,
        double b,
        double c,
        double alpha,
        double beta,
        double gamma,
        const char* spacegroup_name = "P 1"
    );

    [[nodiscard]] static CrystalMetadata from_mosflm_cell(
        const UnitCell& cell,
        const char* spacegroup_name = "P 1"
    );

    [[nodiscard]] UnitCell to_mosflm_cell() const;

    [[nodiscard]] gemmi::UnitCell reciprocal() const {
        return cell.reciprocal();
    }

    [[nodiscard]] std::string describe() const;
};

}  // namespace mosflm::crystal
