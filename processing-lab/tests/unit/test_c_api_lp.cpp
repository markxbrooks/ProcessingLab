#include "mosflm/diffraction/LpCorrection.hpp"
#include "mosflm/legacy/MosflmCApi.h"
#include "test_util.hpp"

using mosflm::DetectorGeometry;

int main() {
    mosflm_detector_c det{};
    det.distance_mm = 100.0;
    det.beam_x_mm = 0.0;
    det.beam_y_mm = 0.0;
    det.pixel_size_x_mm = 0.1;
    det.pixel_size_y_mm = 0.1;
    det.nfast = 2000;
    det.nslow = 2000;
    det.wavelength_A = 1.0;

    DetectorGeometry g;
    g.distance_mm = 100.0;
    g.beam_x_mm = 0.0;
    g.beam_y_mm = 0.0;
    g.pixel_size_x_mm = 0.1;
    g.pixel_size_y_mm = 0.1;
    g.nfast = 2000;
    g.nslow = 2000;
    g.wavelength_A = 1.0;

    // On-axis clamp (no singularity at 2theta = 0).
    CHECK_CLOSE(mosflm_lp_factor(&det, 0.0, 0.0), 1.0, 1e-12);

    // C ABI agrees with the domain implementation off-axis.
    const double xs[] = {50.0, 20.0, -30.0, 0.0};
    const double ys[] = {0.0, 40.0, -10.0, 60.0};
    for (int i = 0; i < 4; ++i) {
        const double c = mosflm_lp_factor(&det, xs[i], ys[i]);
        const double d = mosflm::diffraction::lp_factor(g, xs[i], ys[i]);
        CHECK_CLOSE(c, d, 1e-12);
    }

    // LP is positive and decreases away from the beam (for this geometry).
    const double near_axis = mosflm_lp_factor(&det, 10.0, 0.0);
    const double far_axis = mosflm_lp_factor(&det, 80.0, 0.0);
    CHECK(far_axis > 0.0);
    CHECK(far_axis < near_axis);
    // Radial only: symmetric in x about the beam.
    CHECK_CLOSE(mosflm_lp_factor(&det, -50.0, 0.0),
                mosflm_lp_factor(&det, 50.0, 0.0), 1e-12);

    // apply_lp_inplace scales intensity and sigma by the same factor.
    double intensity = 100.0;
    double sigma = 10.0;
    const double factor = mosflm::diffraction::lp_factor(g, 50.0, 0.0);
    mosflm::diffraction::apply_lp_inplace(g, 50.0, 0.0, intensity, sigma);
    CHECK_CLOSE(intensity, 100.0 * factor, 1e-9);
    CHECK_CLOSE(sigma, 10.0 * factor, 1e-9);

    return mosflm_test::finish("c_api_lp");
}
