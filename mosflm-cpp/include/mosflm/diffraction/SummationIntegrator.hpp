#pragma once

#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/image/DiffractionImage.hpp"

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

struct IntegrateRequest {
    /// Half-width of the peak box in pixels (box is ``2*half+1`` square).
    int box_half_px{3};
    /// Background frame thickness outside the peak box (pixels).
    int bg_frame_px{2};
    /// If true, scale intensity/sigma by unpolarized ``lp_factor``.
    bool apply_lp{false};
};

/// Box summation integrator for predicted spots on a still image.
///
/// Pixel convention: detector origin at a corner; beam centre at
/// ``(nfast/2, nslow/2)`` in pixel space; predicted ``x_mm``/``y_mm`` are
/// beam-relative (as from :class:`ReflectionPredictor`).
class SummationIntegrator {
public:
    [[nodiscard]] std::vector<IntegratedReflection> integrate(
        const image::DiffractionImage& image,
        const DetectorGeometry& detector,
        const std::vector<PredictedReflection>& predictions,
        const IntegrateRequest& request = {}
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
