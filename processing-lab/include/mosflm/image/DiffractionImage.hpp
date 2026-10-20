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

    /// MOSFLM detector omega (``.spt`` header field 5), degrees: spot (x, y)
    /// are rotated by this before the φ rotation about lab Y. 0 = rotation
    /// axis along +slow (ADSC/SMV convention); 270 = axis along +fast.
    /// Derived from NeXus axis/pixel-direction vectors for HDF5 masters.
    double rotation_axis_omega_deg{0.0};

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

/// Dectris Eiger NeXus / ``*_master.h5`` and DLS-style ``.nxs`` reader.
///
/// Requires HDF5 at build time (``MOSFLM_HAS_HDF5``). Registers vendored
/// bitshuffle + LZ4 filters. ``frame`` is 1-based (MOSFLM PROCESS convention).
class Hdf5ImageReader {
public:
    [[nodiscard]] DiffractionImage read(
        const std::filesystem::path& path,
        int frame = 1
    ) const;

    [[nodiscard]] int num_images(const std::filesystem::path& path) const;
};

}  // namespace mosflm::image
