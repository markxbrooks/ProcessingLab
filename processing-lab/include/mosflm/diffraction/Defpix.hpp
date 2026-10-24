#pragma once

#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/image/DiffractionImage.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

namespace mosflm::diffraction {

/// An untrusted detector region in XDS pixel coordinates: 1-based, ``x``
/// along fast (image column + 1), ``y`` along slow (row + 1), bounds
/// inclusive, as ``UNTRUSTED_RECTANGLE=`` etc. in XDS.INP.
struct UntrustedRegion {
    enum class Shape { Rectangle, Ellipse, Quadrilateral };
    Shape shape{Shape::Rectangle};
    /// Rectangle: ``x1 x2 y1 y2``. Ellipse: its bounding box ``x1 x2 y1 y2``.
    /// Quadrilateral: corners ``x1 y1 x2 y2 x3 y3 x4 y4`` in order.
    std::vector<double> values;

    /// Whether pixel (``col``, ``row``), 0-based, lies inside.
    [[nodiscard]] bool contains(int col, int row) const;
    /// ``RECTANGLE``, ``ELLIPSE`` or ``QUADRILATERAL``.
    [[nodiscard]] std::string shape_name() const;
};

/// Parse ``SHAPE v1 v2 …`` (the arguments of CONTROL ``UNTRUSTED``).
///
/// :raises std::invalid_argument: For an unknown shape or the wrong number
///     of values.
[[nodiscard]] UntrustedRegion parse_untrusted_region(const std::vector<std::string>& args);

/// Resolution range (Å) excluded from integration, either order.
struct ResolutionRange {
    double d1_A{0.0};
    double d2_A{0.0};
    [[nodiscard]] bool contains(double d_A) const;
};

/// The ice-ring ranges XDS users commonly exclude (CONTROL ``EXCLUDE ICE``).
[[nodiscard]] std::vector<ResolutionRange> ice_ring_ranges();

/// What DEFPIX decided about a pixel.
enum class PixelStatus : std::uint8_t {
    Trusted = 0,
    /// Masked or overloaded in every sampled frame by the detector itself
    /// (module gaps, the detector's bad-pixel mask).
    DetectorMasked,
    /// Inside an ``UNTRUSTED`` region.
    UntrustedRegion,
    /// In an ``EXCLUDE RESOLUTION`` / ``EXCLUDE ICE`` range.
    ExcludedResolution,
    /// Background outside the trusted value range for its resolution:
    /// shaded by the beamstop or other hardware (XDS ``ABS.cbf`` test).
    Shaded,
    /// Overloaded in at least half the sampled frames (stuck at the counter
    /// maximum), or in two or more and below ``hot_min_counts`` in the rest
    /// (flickering between background and saturation).
    Defective,
    /// High in every sampled frame, far above the background at its
    /// resolution.
    Hot,
};

/// DEFPIX settings (CONTROL ``DEFPIX``, ``UNTRUSTED``, ``EXCLUDE``).
struct DefpixSettings {
    bool enabled{true};
    std::vector<UntrustedRegion> untrusted;
    std::vector<ResolutionRange> excluded_resolution;
    /// Frames sampled across the sweep for the background (at most 50).
    int frames{10};
    /// Detect shaded regions and defective / hot pixels (need ≥ 3 frames).
    bool find_shadows{true};
    bool find_defects{true};
    /// XDS ``VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS``: a block's background
    /// relative to its resolution ring, ×10000, must lie in this range.
    double value_range_low{7000.0};
    double value_range_high{30000.0};
    /// A block is judged only with at least this many counts over the
    /// sampled frames; block size adapts so a typical block has them.
    double shadow_min_counts{400.0};
    /// Hot: every sampled frame at least ``hot_min_counts`` and
    /// ``hot_factor`` × the mean background of its resolution ring.
    double hot_min_counts{10.0};
    double hot_factor{20.0};
};

/// Pixel classification over the sweep.
struct DefpixResult {
    int width{0};
    int height{0};
    /// :class:`PixelStatus` per pixel, row-major.
    std::vector<std::uint8_t> status;
    /// Spot-free background per pixel, counts per frame (the mean over the
    /// sampled frames without each pixel's largest value), row-major; −1
    /// where the pixel never held a count.
    std::vector<float> background;
    std::vector<int> frames_used;
    /// Edge of the square blocks the shadow test judged (px); 0 if not run.
    int shadow_block_px{0};

    [[nodiscard]] std::size_t count(PixelStatus s) const;
    /// ``1`` = use (trusted, or masked by the detector, which the pixel
    /// values already say), ``0`` = rejected by DEFPIX.
    [[nodiscard]] std::vector<std::uint8_t> trusted_mask() const;
    /// Pixels DEFPIX rejected (every status except Trusted / DetectorMasked).
    [[nodiscard]] std::size_t n_rejected() const;
};

/// ``n`` frames spread evenly over ``first``…``last`` (all if fewer).
[[nodiscard]] std::vector<int> defpix_frames(int first, int last, int n);

/// Classify every pixel of a sweep (XDS DEFPIX with INIT's background).
///
/// Frames are read one at a time through ``read_frame``, so memory stays at
/// a few values per pixel. The background is each pixel's sum over the
/// sampled frames without its largest value, so a Bragg spot crossing it
/// in one sampled frame does not count.
///
/// :param read_frame: Returns frame ``n`` (as numbered in ``frames``).
/// :param frames: Frames to sample, e.g. from :func:`defpix_frames`.
/// :param detector: Pixel size, distance and wavelength (for resolution).
/// :param beam_x_mm: Beam position on the panel (mm from the first pixel,
///     fast); with ``beam_y_mm`` the centre of the resolution rings.
/// :param settings: Regions, ranges and detection thresholds.
/// :returns: The classification.
/// :raises std::runtime_error: If the frames differ in size.
[[nodiscard]] DefpixResult run_defpix(
    const std::function<image::DiffractionImage(int)>& read_frame,
    const std::vector<int>& frames,
    const DetectorGeometry& detector,
    double beam_x_mm,
    double beam_y_mm,
    const DefpixSettings& settings
);

/// Write ``result`` as an XDS ``BKGPIX.cbf``-style control image: the
/// background (sum over the sampled frames) for trusted pixels, −3 where
/// DEFPIX rejected the pixel and −1 where the detector masks it.
void write_defpix_cbf(const std::filesystem::path& path, const DefpixResult& result);

}  // namespace mosflm::diffraction
