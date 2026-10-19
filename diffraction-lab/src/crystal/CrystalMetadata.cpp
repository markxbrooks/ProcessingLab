#include "mosflm/crystal/CrystalMetadata.hpp"

#include "mosflm/crystal/SpaceGroupInfo.hpp"

#include <cstdio>

namespace mosflm::crystal {

CrystalMetadata CrystalMetadata::from_parameters(
    double a,
    double b,
    double c,
    double alpha,
    double beta,
    double gamma,
    const char* spacegroup_name
) {
    CrystalMetadata meta;
    meta.cell = gemmi::UnitCell(a, b, c, alpha, beta, gamma);
    if (!meta.cell.is_crystal()) {
        throw std::invalid_argument("CrystalMetadata: invalid unit cell");
    }
    meta.spacegroup = SpaceGroupInfo::find_by_name(spacegroup_name);
    if (meta.spacegroup == nullptr) {
        throw std::invalid_argument(
            std::string("Unknown space group: ") + spacegroup_name
        );
    }
    return meta;
}

CrystalMetadata CrystalMetadata::from_mosflm_cell(
    const UnitCell& cell,
    const char* spacegroup_name
) {
    return from_parameters(
        cell.a,
        cell.b,
        cell.c,
        cell.alpha_deg,
        cell.beta_deg,
        cell.gamma_deg,
        spacegroup_name
    );
}

UnitCell CrystalMetadata::to_mosflm_cell() const {
    UnitCell out;
    out.a = cell.a;
    out.b = cell.b;
    out.c = cell.c;
    out.alpha_deg = cell.alpha;
    out.beta_deg = cell.beta;
    out.gamma_deg = cell.gamma;
    return out;
}

std::string CrystalMetadata::describe() const {
    char buf[192];
    const std::string sg =
        (spacegroup != nullptr) ? spacegroup->xhm() : std::string("(none)");
    std::snprintf(
        buf,
        sizeof(buf),
        "a=%.4f b=%.4f c=%.4f α=%.2f β=%.2f γ=%.2f sg=%s",
        cell.a,
        cell.b,
        cell.c,
        cell.alpha,
        cell.beta,
        cell.gamma,
        sg.c_str()
    );
    return std::string(buf);
}

}  // namespace mosflm::crystal
