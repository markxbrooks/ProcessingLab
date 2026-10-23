#pragma once

#include "mosflm/io/ImportedSweepGeometry.hpp"

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace mosflm::io {

/// Parse Xia2 ``.xinfo``: one sweep (the first, or a named one) and the
/// wavelength it refers to (else the first ``WAVELENGTH`` block).
class XinfoImporter {
public:
    [[nodiscard]] static ImportedSweepGeometry parse_string(
        std::string_view text
    );
    /// The sweep named ``sweep`` (``BEGIN SWEEP <name>``).
    /// :raises std::runtime_error: If there is no such sweep.
    [[nodiscard]] static ImportedSweepGeometry parse_string(
        std::string_view text, std::string_view sweep
    );
    [[nodiscard]] static ImportedSweepGeometry parse_file(
        const std::filesystem::path& path
    );
    [[nodiscard]] static ImportedSweepGeometry parse_file(
        const std::filesystem::path& path, std::string_view sweep
    );
    /// Sweep names in file order.
    [[nodiscard]] static std::vector<std::string> sweep_names(std::string_view text);
    [[nodiscard]] static std::vector<std::string> sweep_names_in_file(
        const std::filesystem::path& path
    );
};

}  // namespace mosflm::io
