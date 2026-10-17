#pragma once

#include "mosflm/crystal/CrystalMetadata.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace mosflm::crystal {

struct SimpleReflection {
    int h{};
    int k{};
    int l{};
    float intensity{};
    float sigma{};
};

/// Minimal MTZ writer/reader for integrate comparisons (Gemmi Mtz API).
class MtzAdapter {
public:
    /// Write a single-dataset MTZ with H K L I SIGI columns.
    static void write_simple(
        const std::filesystem::path& path,
        const CrystalMetadata& crystal,
        const std::vector<SimpleReflection>& reflections
    );

    /// Read H K L I SIGI from an MTZ written by :func:`write_simple`.
    [[nodiscard]] static std::vector<SimpleReflection> read_simple(
        const std::filesystem::path& path
    );
};

}  // namespace mosflm::crystal
