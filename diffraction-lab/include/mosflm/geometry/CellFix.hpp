#pragma once

#include "mosflm/geometry/UnitCell.hpp"

#include <gemmi/symmetry.hpp>

#include <string>

namespace mosflm {

/// MOSFLM ``ICRYST`` crystal-system codes (CELLFIX / CELLCHK).
enum class CrystalSystemCode {
    Triclinic = 1,
    Monoclinic = 2,
    Orthorhombic = 3,
    Tetragonal = 4,
    Trigonal = 5,
    Hexagonal = 6,
    Cubic = 7,
    Rhombohedral = 8
};

/// Map a gemmi space group to MOSFLM ``ICRYST`` (1–8).
[[nodiscard]] CrystalSystemCode icryst_from_spacegroup(
    const gemmi::SpaceGroup* sg
);

/// MOSFLM ``CELLFIX``: project cell edges/angles onto the crystal-system
/// metric. ``softness`` in (0,1] blends toward the ideal (1 = hard fix).
[[nodiscard]] UnitCell cellfix(
    const UnitCell& cell,
    CrystalSystemCode icryst,
    double softness = 1.0
);

/// Convenience: CELLFIX using the space-group crystal system.
[[nodiscard]] UnitCell cellfix_for_spacegroup(
    const UnitCell& cell,
    const gemmi::SpaceGroup* sg,
    double softness = 1.0
);

/// Result of comparing the cell metric to the assigned space group.
struct LatticeConsistency {
    bool ok{true};
    std::string message;
    /// Laue-group HM from gemmi lattice symmetry (empty if unavailable).
    std::string lattice_laue_hm;
    int lattice_nops{0};
    int spacegroup_nops{0};
};

/// Use gemmi ``find_lattice_symmetry`` to check whether the cell metric
/// supports at least the Laue symmetry implied by ``sg``.
[[nodiscard]] LatticeConsistency check_lattice_vs_spacegroup(
    const UnitCell& cell,
    const gemmi::SpaceGroup* sg,
    double max_obliq_deg = 2.0
);

}  // namespace mosflm
