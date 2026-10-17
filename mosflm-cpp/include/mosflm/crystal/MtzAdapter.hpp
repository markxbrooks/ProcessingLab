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

/// Minimal MTZ writer for future integrate comparisons (Gemmi Mtz API).
class MtzAdapter {
public:
    /// Write a single-dataset MTZ with H K L I SIGI columns.
    static void write_simple(
        const std::filesystem::path& path,
        const CrystalMetadata& crystal,
        const std::vector<SimpleReflection>& reflections
    );
};

}  // namespace mosflm::crystal
