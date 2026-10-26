#pragma once

#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/geometry/LabFrame.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace mosflm::image {

/// Pixel values below 0 are not measurements. Every reader maps its
/// detector's own markers onto these two; integration skips both.
/// Masked: module gap, bad or untrusted pixel (Pilatus CBF −1 and −2, Eiger
/// 2ⁿ−1). Overloaded: at or above ``ImageMetadata::saturation_counts``.
inline constexpr std::int32_t kPixelMasked = -1;
inline constexpr std::int32_t kPixelOverloaded = -2;

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

    /// Overload threshold (counts) from the header; 0 if unknown.
    double saturation_counts{0.0};

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
        // Detector omega → MOSFLM CCOM (INVERTX convention, i.e. a MOSFLM
        // A-matrix). CONTROL overrides this when AUTOINDEX used a
        // non-inverted ProcessingLab ``.spt`` (see ControlState::spt_invertx).
        g.ccom_deg = ccom_from_spt_omega(rotation_axis_omega_deg, true);
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

/// Write ``image`` as a miniCBF (byte-offset compressed, as Pilatus /
/// Eiger CBFs are), readable by :class:`CbfImageReader`, XDS-Viewer, adxv
/// and fabio.
///
/// :param path: Destination file.
/// :param image: Pixels (row-major, ``width`` × ``height``) and size.
/// :param header_contents: ``_array_data.header_contents`` text (``# Key
///     value`` lines), or empty.
/// :raises std::runtime_error: If the file cannot be written.
void write_minicbf(
    const std::filesystem::path& path,
    const DiffractionImage& image,
    const std::string& header_contents = {}
);

/// ADSC / SMV reader (ASCII ``{KEY=VAL;}`` header + unsigned_short pixels).
///
/// Enough for the CCP4 insulin example (``insulin_1_001.img``). Supports
/// little-endian ``TYPE=unsigned_short`` with ``HEADER_BYTES`` / ``SIZE1`` /
/// ``SIZE2``. Pixels are stored transposed (``width = SIZE2``): the fast
/// index runs along the file's slow axis, so (col, row) = MOSFLM
/// (XDET, YDET) and the rotation axis (file fast axis) is along +slow.
class SmvImageReader {
public:
    [[nodiscard]] DiffractionImage read(const std::filesystem::path& path) const;
};

/// Dectris Eiger NeXus / ``*_master.h5`` and DLS-style ``.nxs`` reader.
///
/// Requires HDF5 at build time (``MOSFLM_HAS_HDF5``). Registers vendored
/// bitshuffle + LZ4 filters, and decodes one-frame bitshuffle-LZ4 chunks
/// itself. Reads DECTRIS ``data_NNNNNN`` blocks or one NXmx ``data`` dataset,
/// applies ``pixel_mask`` (−1) and honours ``units`` (adapted in part from
/// durin; see third_party/durin/LICENSE). ``frame`` is 1-based (MOSFLM PROCESS convention).
class Hdf5ImageReader {
public:
    [[nodiscard]] DiffractionImage read(
        const std::filesystem::path& path,
        int frame = 1
    ) const;

    [[nodiscard]] int num_images(const std::filesystem::path& path) const;
};

}  // namespace mosflm::image
