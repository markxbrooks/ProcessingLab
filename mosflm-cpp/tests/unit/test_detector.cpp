#include "mosflm/geometry/DetectorGeometry.hpp"
#include "test_util.hpp"

using mosflm::DetectorGeometry;

int main() {
    DetectorGeometry det;
    det.distance_mm = 100.0;
    det.beam_x_mm = 0.0;
    det.beam_y_mm = 0.0;
    det.pixel_size_x_mm = 0.1;
    det.pixel_size_y_mm = 0.1;
    det.nfast = 2000;
    det.nslow = 2000;
    det.wavelength_A = 1.0;

    // Ray exactly along the beam axis hits the beam centre.
    double x = 0.0, y = 0.0;
    CHECK((det.lab_to_detector_mm({0.0, 0.0, 1.0}, x, y)));
    CHECK_CLOSE(x, 0.0, 1e-12);
    CHECK_CLOSE(y, 0.0, 1e-12);

    // 45 degrees off axis lands at x = distance.
    CHECK((det.lab_to_detector_mm({0.5, 0.0, 0.5}, x, y)));
    CHECK_CLOSE(x, 100.0, 1e-9);
    CHECK_CLOSE(y, 0.0, 1e-9);

    // Back-scattered ray (negative z) cannot hit the forward detector.
    CHECK((!det.lab_to_detector_mm({0.0, 0.0, -1.0}, x, y)));

    // Pixel mapping.
    double xf = 0.0, yf = 0.0;
    det.mm_to_pixel(10.0, -5.0, xf, yf);
    CHECK_CLOSE(xf, 100.0, 1e-12);
    CHECK_CLOSE(yf, -50.0, 1e-12);

    // two_theta: on-axis -> 0; radial 100 mm at 100 mm distance -> 45 deg.
    CHECK_CLOSE(det.two_theta_deg(0.0, 0.0), 0.0, 1e-9);
    CHECK_CLOSE(det.two_theta_deg(100.0, 0.0), 45.0, 1e-9);
    CHECK_CLOSE(det.two_theta_deg(0.0, 100.0), 45.0, 1e-9);

    // Off-centre beam shifts the origin.
    det.beam_x_mm = 10.0;
    CHECK((det.lab_to_detector_mm({0.0, 0.0, 1.0}, x, y)));
    CHECK_CLOSE(x, 10.0, 1e-12);

    // YSCALE adjusts YG ↔ C++ x_mm pixel pitch in mm_to_pixel.
    det.beam_x_mm = 0.0;
    det.yscale = 2.0;
    det.mm_to_pixel(10.0, 10.0, xf, yf);
    CHECK_CLOSE(xf, 200.0, 1e-12);
    CHECK_CLOSE(yf, 100.0, 1e-12);

    // First-order TWIST couples YG into XG (→ C++ y_mm).
    det.yscale = 1.0;
    det.twist_mrad = 1000.0;  // 1 rad
    CHECK((det.lab_to_detector_mm({0.5, 0.0, 0.5}, x, y)));
    CHECK_CLOSE(x, 100.0, 1e-9);
    CHECK_CLOSE(y, 100.0, 1e-9);

    return mosflm_test::finish("detector");
}
