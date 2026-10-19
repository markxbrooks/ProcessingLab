#include "mosflm/geometry/CellFix.hpp"
#include "test_util.hpp"

#include <gemmi/symmetry.hpp>

#include <cmath>
#include <iostream>

int main() {
    using mosflm::UnitCell;
    using mosflm::cellfix_for_spacegroup;
    using mosflm::check_lattice_vs_spacegroup;

    const auto* p4 = gemmi::find_spacegroup_by_number(75);
    CHECK(p4 != nullptr);

    const UnitCell observed{57.98, 58.04, 148.16, 90.00, 90.30, 90.10};
    const UnitCell fixed = cellfix_for_spacegroup(observed, p4, 1.0);
    CHECK_MSG(std::fabs(fixed.a - fixed.b) < 1e-9, "tetragonal a!=b after CELLFIX");
    CHECK_MSG(std::fabs(fixed.a - 58.01) < 1e-6, "a,b mean");
    CHECK_MSG(std::fabs(fixed.alpha_deg - 90.0) < 1e-9, "alpha");
    CHECK_MSG(std::fabs(fixed.beta_deg - 90.0) < 1e-9, "beta");
    CHECK_MSG(std::fabs(fixed.gamma_deg - 90.0) < 1e-9, "gamma");
    CHECK_MSG(std::fabs(fixed.c - 148.16) < 1e-9, "c unchanged");

    const auto chk = check_lattice_vs_spacegroup(observed, p4, 2.0);
    // Near-tetragonal metric should still find lattice ops; warn on |a-b|.
    CHECK_MSG(!chk.ok, "expected warn for unconstrained tetragonal cell");
    std::cout << "cell_fix: " << chk.message << '\n';

    const UnitCell already{58.01, 58.01, 148.16, 90, 90, 90};
    const auto chk2 = check_lattice_vs_spacegroup(already, p4, 2.0);
    CHECK_MSG(chk2.ok, "ideal tetragonal should pass");

    return mosflm_test::finish("cell_fix");
}
