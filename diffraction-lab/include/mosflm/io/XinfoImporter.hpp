#pragma once

#include "mosflm/io/ImportedSweepGeometry.hpp"

#include <filesystem>
#include <string>
#include <string_view>

namespace mosflm::io {

/// Parse Xia2 ``.xinfo`` (first wavelength + first sweep only).
class XinfoImporter {
public:
    [[nodiscard]] static ImportedSweepGeometry parse_string(
        std::string_view text
    );

    [[nodiscard]] static ImportedSweepGeometry parse_file(
        const std::filesystem::path& path
    );
};

}  // namespace mosflm::io
