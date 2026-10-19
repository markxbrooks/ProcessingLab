#pragma once

#include "mosflm/io/ImportedSweepGeometry.hpp"

#include <filesystem>
#include <string>
#include <string_view>

namespace mosflm::io {

/// Parse an XDS ``XDS.INP``-style ``KEY=value`` file into CONTROL geometry.
///
/// Does not read orientation (use GXPARM/XPARM separately). Omits ``BEAM`` when
/// ``ORGX``/``ORGY`` are absent.
class XdsInpImporter {
public:
    [[nodiscard]] static ImportedSweepGeometry parse_string(
        std::string_view text
    );

    [[nodiscard]] static ImportedSweepGeometry parse_file(
        const std::filesystem::path& path
    );
};

}  // namespace mosflm::io
