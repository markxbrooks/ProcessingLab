#include "mosflm/geometry/CellFix.hpp"

#include <gemmi/twin.hpp>

#include <algorithm>
#include <cmath>
#include <sstream>

namespace mosflm {
namespace {

double soften(double value, double target, double softness) {
    const double s = std::clamp(softness, 0.0, 1.0);
    return value - (value - target) * s;
}

}  // namespace

CrystalSystemCode icryst_from_spacegroup(const gemmi::SpaceGroup* sg) {
    if (sg == nullptr) {
        return CrystalSystemCode::Triclinic;
    }
    switch (sg->crystal_system()) {
    case gemmi::CrystalSystem::Triclinic:
        return CrystalSystemCode::Triclinic;
    case gemmi::CrystalSystem::Monoclinic:
        return CrystalSystemCode::Monoclinic;
    case gemmi::CrystalSystem::Orthorhombic:
        return CrystalSystemCode::Orthorhombic;
    case gemmi::CrystalSystem::Tetragonal:
        return CrystalSystemCode::Tetragonal;
    case gemmi::CrystalSystem::Trigonal: {
        // Hall starting with R → rhombohedral axes (ICRYST 8).
        const std::string hall(sg->hall);
        if (!hall.empty() && (hall[0] == 'R' || hall[0] == 'r')) {
            return CrystalSystemCode::Rhombohedral;
        }
        return CrystalSystemCode::Trigonal;
    }
    case gemmi::CrystalSystem::Hexagonal:
        return CrystalSystemCode::Hexagonal;
    case gemmi::CrystalSystem::Cubic:
        return CrystalSystemCode::Cubic;
    }
    return CrystalSystemCode::Triclinic;
}

UnitCell cellfix(
    const UnitCell& cell,
    CrystalSystemCode icryst,
    double softness
) {
    UnitCell out = cell;
    const double s = std::clamp(softness, 0.0, 1.0);
    if (s <= 0.0 || icryst == CrystalSystemCode::Triclinic) {
        return out;
    }

    switch (icryst) {
    case CrystalSystemCode::Monoclinic:
        // Soften alpha and gamma toward 90° (beta free).
        out.alpha_deg = soften(out.alpha_deg, 90.0, s);
        out.gamma_deg = soften(out.gamma_deg, 90.0, s);
        break;
    case CrystalSystemCode::Orthorhombic:
        out.alpha_deg = soften(out.alpha_deg, 90.0, s);
        out.beta_deg = soften(out.beta_deg, 90.0, s);
        out.gamma_deg = soften(out.gamma_deg, 90.0, s);
        break;
    case CrystalSystemCode::Tetragonal: {
        out.alpha_deg = soften(out.alpha_deg, 90.0, s);
        out.beta_deg = soften(out.beta_deg, 90.0, s);
        out.gamma_deg = soften(out.gamma_deg, 90.0, s);
        const double ab = 0.5 * (out.a + out.b);
        out.a = out.a + (ab - out.a) * s;
        out.b = out.b + (ab - out.b) * s;
        break;
    }
    case CrystalSystemCode::Trigonal:
    case CrystalSystemCode::Hexagonal: {
        out.alpha_deg = soften(out.alpha_deg, 90.0, s);
        out.beta_deg = soften(out.beta_deg, 90.0, s);
        const double d60 = std::fabs(out.gamma_deg - 60.0);
        const double d120 = std::fabs(out.gamma_deg - 120.0);
        const double gtarget = (d60 < d120) ? 60.0 : 120.0;
        out.gamma_deg = soften(out.gamma_deg, gtarget, s);
        const double ab = 0.5 * (out.a + out.b);
        out.a = out.a + (ab - out.a) * s;
        out.b = out.b + (ab - out.b) * s;
        break;
    }
    case CrystalSystemCode::Cubic: {
        out.alpha_deg = soften(out.alpha_deg, 90.0, s);
        out.beta_deg = soften(out.beta_deg, 90.0, s);
        out.gamma_deg = soften(out.gamma_deg, 90.0, s);
        const double m = (out.a + out.b + out.c) / 3.0;
        out.a = out.a + (m - out.a) * s;
        out.b = out.b + (m - out.b) * s;
        out.c = out.c + (m - out.c) * s;
        break;
    }
    case CrystalSystemCode::Rhombohedral: {
        const double m = (out.a + out.b + out.c) / 3.0;
        out.a = out.a + (m - out.a) * s;
        out.b = out.b + (m - out.b) * s;
        out.c = out.c + (m - out.c) * s;
        const double ang = (out.alpha_deg + out.beta_deg + out.gamma_deg) / 3.0;
        out.alpha_deg = out.alpha_deg + (ang - out.alpha_deg) * s;
        out.beta_deg = out.beta_deg + (ang - out.beta_deg) * s;
        out.gamma_deg = out.gamma_deg + (ang - out.gamma_deg) * s;
        break;
    }
    case CrystalSystemCode::Triclinic:
        break;
    }
    return out;
}

UnitCell cellfix_for_spacegroup(
    const UnitCell& cell,
    const gemmi::SpaceGroup* sg,
    double softness
) {
    return cellfix(cell, icryst_from_spacegroup(sg), softness);
}

LatticeConsistency check_lattice_vs_spacegroup(
    const UnitCell& cell,
    const gemmi::SpaceGroup* sg,
    double max_obliq_deg
) {
    LatticeConsistency out;
    if (sg == nullptr) {
        out.ok = true;
        out.message = "no space group";
        return out;
    }
    gemmi::UnitCell gcell(
        cell.a, cell.b, cell.c, cell.alpha_deg, cell.beta_deg, cell.gamma_deg
    );
    if (!gcell.is_crystal()) {
        out.ok = false;
        out.message = "invalid unit cell";
        return out;
    }
    char centring = sg->centring_type();
    if (centring == '\0') {
        centring = 'P';
    }
    const gemmi::GroupOps lops =
        gemmi::find_lattice_symmetry(gcell, centring, max_obliq_deg);
    out.lattice_nops = static_cast<int>(lops.sym_ops.size());
    out.spacegroup_nops = sg->operations().order();
    const gemmi::SpaceGroup* lsg = gemmi::find_spacegroup_by_ops(lops);
    if (lsg != nullptr) {
        out.lattice_laue_hm = lsg->hm;
    }
    // Laue group of assigned SG should not exceed metric lattice ops.
    // Compare crystal systems: warn if SG needs higher metric than cell.
    const auto sg_sys = sg->crystal_system();
    bool metric_ok = true;
    std::ostringstream msg;
    if (sg_sys == gemmi::CrystalSystem::Tetragonal ||
        sg_sys == gemmi::CrystalSystem::Hexagonal ||
        sg_sys == gemmi::CrystalSystem::Cubic) {
        if (out.lattice_nops < 4) {
            metric_ok = false;
            msg << "cell metric (lattice nops=" << out.lattice_nops
                << ") weaker than " << sg->hm;
        }
    }
    if (sg_sys == gemmi::CrystalSystem::Tetragonal) {
        const double dab = std::fabs(cell.a - cell.b);
        const double dang = std::max(
            {std::fabs(cell.alpha_deg - 90.0), std::fabs(cell.beta_deg - 90.0),
             std::fabs(cell.gamma_deg - 90.0)}
        );
        // Tight enough that unconstrained autoindex cells (Δab~0.05 Å,
        // Δangle~0.3°) still prompt CELLFIX before PROCESS/MTZ.
        if (dab > 0.02 || dang > 0.05) {
            metric_ok = false;
            if (!msg.str().empty()) {
                msg << "; ";
            }
            msg << "|a-b|=" << dab << " max|angle-90|=" << dang
                << " (tetragonal expects a=b, angles≈90)";
        }
    }
    out.ok = metric_ok;
    if (msg.str().empty()) {
        msg << "lattice " << (lsg ? lsg->hm : "?") << " nops=" << out.lattice_nops
            << " vs SG " << sg->hm << " nops=" << out.spacegroup_nops;
    }
    out.message = msg.str();
    return out;
}

}  // namespace mosflm
