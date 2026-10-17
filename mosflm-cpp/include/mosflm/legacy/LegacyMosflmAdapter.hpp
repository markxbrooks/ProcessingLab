#pragma once

#include "mosflm/geometry/UnitCell.hpp"

#include <string>

namespace mosflm::legacy {

/// Narrow façade for future Fortran / C interoperability.
///
/// Keep COMMON blocks, column-major arrays, and compiler mangling **inside**
/// this layer. Domain code should never see raw ``CELL(6)`` buffers.
class LegacyMosflmAdapter {
public:
    /// Pack a :class:`UnitCell` into MOSFLM ``CELL(6)`` order.
    static void to_cell6(const UnitCell& cell, float out[6]);

    /// Unpack MOSFLM ``CELL(6)`` into a :class:`UnitCell`.
    static UnitCell from_cell6(const float in[6]);

    /// Human-readable note about the linked oracle binary, if any.
    [[nodiscard]] static std::string oracle_hint();
};

}  // namespace mosflm::legacy
