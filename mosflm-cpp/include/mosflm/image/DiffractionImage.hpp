#pragma once

#include "mosflm/geometry/DetectorGeometry.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace mosflm::image {

/// Experiment metadata normalized for MOSFLM-style processing.
struct ImageMetadata {
    int width{};
    int height{};

    double pixel_size_x_mm{};
    double pixel_size_y_mm{};

    double wavelength_angstrom{};
    double detector_distance_mm{};

    double beam_x_mm{};
    double beam_y_mm{};

    double oscillation_start_deg{};
    double oscillation_range_deg{};

    std::string detector_id;

    /// Map into the geometry layer used by prediction / integration.
    [[nodiscard]] DetectorGeometry to_detector_geometry() const {
        DetectorGeometry g;
        g.nfast = width;
        g.nslow = height;
        g.pixel_size_x_mm = pixel_size_x_mm;
        g.pixel_size_y_mm = pixel_size_y_mm;
        g.distance_mm = detector_distance_mm;
        g.beam_x_mm = beam_x_mm;
        g.beam_y_mm = beam_y_mm;
        g.wavelength_A = wavelength_angstrom;
        return g;
    }
};

struct DiffractionImage {
    ImageMetadata metadata;
    std::vector<std::int32_t> pixels;
};

/// CBFlib-backed reader (Option A). See docs/CBF_STATUS.md.
class CbfImageReader {
public:
    [[nodiscard]] DiffractionImage read(const std::filesystem::path& path) const;
};

/// ADSC / SMV reader (ASCII ``{KEY=VAL;}`` header + unsigned_short pixels).
///
/// Enough for the CCP4 insulin example (``insulin_1_001.img``). Supports
/// little-endian ``TYPE=unsigned_short`` with ``HEADER_BYTES`` / ``SIZE1`` /
/// ``SIZE2``.
class SmvImageReader {
public:
    [[nodiscard]] DiffractionImage read(const std::filesystem::path& path) const;
};

}  // namespace mosflm::image
