#include "mosflm/crystal/SpaceGroupInfo.hpp"

namespace mosflm::crystal {

const gemmi::SpaceGroup* SpaceGroupInfo::find_by_name(const char* name) {
    return gemmi::find_spacegroup_by_name(name);
}

std::string SpaceGroupInfo::hall(const gemmi::SpaceGroup* sg) {
    if (sg == nullptr) {
        return {};
    }
    return std::string(sg->hall);
}

int SpaceGroupInfo::number(const gemmi::SpaceGroup* sg) {
    if (sg == nullptr) {
        return 0;
    }
    return sg->number;
}

}  // namespace mosflm::crystal
