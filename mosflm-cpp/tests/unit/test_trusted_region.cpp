#include "mosflm/diffraction/TrustedRegion.hpp"
#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/image/DiffractionImage.hpp"
#include "test_util.hpp"

#include <iostream>

int main() {
    using mosflm::DetectorGeometry;
    using mosflm::diffraction::TrustedRegionSettings;
    using mosflm::diffraction::build_trusted_mask;
    using mosflm::diffraction::is_trusted;
    using mosflm::image::DiffractionImage;

    DiffractionImage img;
    img.metadata.width = 100;
    img.metadata.height = 80;
    img.pixels.assign(static_cast<std::size_t>(100 * 80), 0);

    DetectorGeometry det;
    det.pixel_size_x_mm = 0.075;
    det.pixel_size_y_mm = 0.075;
    det.distance_mm = 200.0;
    det.wavelength_A = 1.0;
    det.beam_x_mm = 3.75;  // col 50
    det.beam_y_mm = 3.0;   // row 40

    TrustedRegionSettings s;
    s.trusted_radius_frac = 0.5;
    s.beam_x_mm = det.beam_x_mm;
    s.beam_y_mm = det.beam_y_mm;
    const auto mask = build_trusted_mask(img, det, s);
    CHECK_EQ(static_cast<int>(mask.size()), 8000);
    CHECK(is_trusted(mask, 100, 50, 40));
    CHECK(!is_trusted(mask, 100, 0, 0));
    CHECK(!is_trusted(mask, 100, 99, 79));

    return mosflm_test::finish("trusted_region");
}
