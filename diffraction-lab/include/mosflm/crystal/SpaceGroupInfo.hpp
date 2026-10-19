#pragma once

#include <gemmi/symmetry.hpp>

#include <string>

namespace mosflm::crystal {

/// Thin helpers around Gemmi space-group tables.
struct SpaceGroupInfo {
    [[nodiscard]] static const gemmi::SpaceGroup* find_by_name(const char* name);

    [[nodiscard]] static const gemmi::SpaceGroup* find_by_number(int number);

    [[nodiscard]] static std::string hall(const gemmi::SpaceGroup* sg);

    [[nodiscard]] static int number(const gemmi::SpaceGroup* sg);
};

}  // namespace mosflm::crystal
