#pragma once

#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/image/DiffractionImage.hpp"

#include <span>
#include <vector>

namespace mosflm::diffraction {

/// One summation-integrated observation (no profile fitting).
struct IntegratedReflection {
    int h{};
    int k{};
    int l{};
    double x_mm{};
    double y_mm{};
    double intensity{};
    double sigma{};
    int npix_peak{};
    int npix_bg{};
    bool ok{false};
};

/// Algorithm options for box summation integrate.
struct IntegrationSettings {
    /// Half-width of the peak box in pixels (box is ``2*half+1`` square).
    int box_half_px{3};
    /// Background frame thickness outside the peak box (pixels).
    int bg_frame_px{2};
    /// If true, scale intensity/sigma by ``lp_factor``.
    bool apply_lp{false};
    /// XDS ``FRACTION_OF_POLARIZATION`` (0 = unpolarized; ~0.98 synchrotron).
    double polarization_fraction{0.0};
    /// If true, map beam-relative mm to ADSC/MOSFLM image pixels via
    /// :func:`cpp_mm_to_mosflm_xdet` (``XG↔y_mm``, ``YG↔x_mm``) using
    /// ``image_beam_*_mm``. Default false keeps the synthetic centre-origin
    /// convention used by vertical-slice tests.
    bool mosflm_detector_axes{false};
    double image_beam_x_mm{0.0};
    double image_beam_y_mm{0.0};
};

/// Inputs for one summation-integrate call.
///
/// ``image``, ``detector``, and ``reflections`` are referenced (not copied);
/// they must outlive the call.
struct IntegrationRequest {
    const image::DiffractionImage& image;
    const DetectorGeometry& detector;
    std::span<const PredictedReflection> reflections;
    IntegrationSettings settings{};
};

/// Box summation integrator for predicted spots on a still image.
///
/// Pixel convention: detector origin at a corner; beam centre at
/// ``(nfast/2, nslow/2)`` in pixel space; predicted ``x_mm``/``y_mm`` are
/// beam-relative (as from :class:`ReflectionPredictor`).
class SummationIntegrator {
public:
    [[nodiscard]] std::vector<IntegratedReflection> integrate(
        const IntegrationRequest& request
    ) const;

    /// Map beam-relative detector mm → integer pixel (column=fast, row=slow).
    static void mm_to_image_pixel(
        const DetectorGeometry& detector,
        double x_mm,
        double y_mm,
        int& col,
        int& row
    );
};

}  // namespace mosflm::diffraction
