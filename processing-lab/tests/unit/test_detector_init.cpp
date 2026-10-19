#include "mosflm/diffraction/DetectorInit.hpp"
#include "mosflm/diffraction/LpCorrection.hpp"
#include "test_util.hpp"

#include <cmath>
#include <iostream>
#include <vector>

int main() {
    using mosflm::DetectorGeometry;
    using mosflm::diffraction::estimate_blank_gain;
    using mosflm::diffraction::lp_factor;
    using mosflm::image::DiffractionImage;

    DiffractionImage img;
    img.metadata.width = 64;
    img.metadata.height = 64;
    img.pixels.assign(64 * 64, 100);
    // Hot corner to ensure samples exist.
    for (int r = 0; r < 8; ++r) {
        for (int c = 0; c < 8; ++c) {
            img.pixels[static_cast<std::size_t>(r) * 64U +
                       static_cast<std::size_t>(c)] = 120;
        }
    }
    const auto est = estimate_blank_gain(img, 8);
    CHECK(est.ok);
    CHECK_MSG(est.blank_adu > 50.0, "blank");
    CHECK_MSG(est.gain > 0.0, "gain");

    DetectorGeometry det;
    det.distance_mm = 100.0;
    det.beam_x_mm = 0.0;
    det.beam_y_mm = 0.0;
    const double unpol = lp_factor(det, 50.0, 0.0, 0.0);
    const double pol = lp_factor(det, 50.0, 0.0, 0.98);
    CHECK_MSG(std::fabs(unpol - pol) > 1e-6, "polarized LP should differ");
    // Horizontal spot: ψ=0 → cos(2ψ)=1 → P smaller than unpolarized.
    CHECK_MSG(pol < unpol, "horizontal pol reduces P");

    return mosflm_test::finish("detector_init");
}
