#include "mosflm/diffraction/LpCorrection.hpp"

#include "mosflm/legacy/MosflmCApi.h"

#include <cmath>

namespace mosflm::diffraction {

double lp_factor(
    const DetectorGeometry& detector,
    double x_mm,
    double y_mm
) {
    mosflm_detector_c d{};
    d.distance_mm = detector.distance_mm;
    d.beam_x_mm = detector.beam_x_mm;
    d.beam_y_mm = detector.beam_y_mm;
    d.pixel_size_x_mm = detector.pixel_size_x_mm;
    d.pixel_size_y_mm = detector.pixel_size_y_mm;
    d.nfast = detector.nfast;
    d.nslow = detector.nslow;
    d.wavelength_A = detector.wavelength_A;
    return mosflm_lp_factor(&d, x_mm, y_mm);
}

void apply_lp_inplace(
    const DetectorGeometry& detector,
    double x_mm,
    double y_mm,
    double& intensity,
    double& sigma
) {
    const double lp = lp_factor(detector, x_mm, y_mm);
    intensity *= lp;
    sigma *= lp;
}

}  // namespace mosflm::diffraction
