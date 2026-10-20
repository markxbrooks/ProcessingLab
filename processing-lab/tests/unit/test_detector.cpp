#include "mosflm/geometry/DetectorGeometry.hpp"
#include "test_util.hpp"

#include <cmath>

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

    // Old MMTOPX: TWIST scales YG via DXTOF = XTOFRA + TWIST*Y (Y in 10 µm).
    // Ray → xg=0, yg=100; ITWIST=1000 → DXTOF = 1 + 1000*DTOR/(10000*100)*10000.
    det.yscale = 1.0;
    det.twist_mrad = 1000.0;
    CHECK((det.lab_to_detector_mm({0.5, 0.0, 0.5}, x, y)));
    {
        constexpr double kDtor = 0.0174532925;
        const double fdist = kDtor / (100.0 * 100.0 * 100.0);
        const double dxtof = 1.0 + 1000.0 * fdist * 10000.0;
        CHECK_CLOSE(x, dxtof * 100.0, 1e-9);
        CHECK_CLOSE(y, 0.0, 1e-9);
    }

    // ROTMAT2 flat normal.
    {
        const Eigen::Vector3d n = mosflm::rotmat2_detnor(0.0, 0.0);
        CHECK_CLOSE(n.x(), 0.0, 1e-15);
        CHECK_CLOSE(n.y(), 0.0, 1e-15);
        CHECK_CLOSE(n.z(), 1.0, 1e-15);
    }

    // NEWT path with zero tilt matches flat (xtofra=1).
    det.twist_mrad = 0.0;
    det.use_rotmat2 = true;
    CHECK((det.lab_to_detector_mm({0.5, 0.0, 0.5}, x, y)));
    CHECK_CLOSE(x, 100.0, 1e-9);
    CHECK_CLOSE(y, 0.0, 1e-9);

    // XTOFRA scales ideal mm.
    det.use_rotmat2 = false;
    det.xtofra = 1.1;
    CHECK((det.lab_to_detector_mm({0.5, 0.0, 0.5}, x, y)));
    CHECK_CLOSE(x, 110.0, 1e-9);
    CHECK_CLOSE(y, 0.0, 1e-9);

    // Farthest Eiger corner is ~1.15 Å; a 24 px inset is slightly lower.
    {
        DetectorGeometry eiger;
        eiger.distance_mm = 199.78169135;
        eiger.wavelength_A = 0.97949913;
        eiger.pixel_size_x_mm = 0.075;
        eiger.pixel_size_y_mm = 0.075;
        eiger.nfast = 4148;
        eiger.nslow = 4362;
        const double beam_x = 166.15771234;
        const double beam_y = 172.97150533;
        const double d_edge =
            mosflm::detector_corner_dmin_A(eiger, beam_x, beam_y, 0.0);
        const double d_inset =
            mosflm::detector_corner_dmin_A(eiger, beam_x, beam_y, 24.0);
        CHECK(d_edge > 1.14 && d_edge < 1.17);
        CHECK(d_inset > d_edge);
        CHECK(d_inset > 1.15 && d_inset < 1.20);

        // Beam-relative mm of panel pixel (0, 0): XDET uses y_mm, YDET uses x_mm.
        CHECK(mosflm::spot_on_panel(
            eiger, -beam_y, -beam_x, 0.0, true, beam_x, beam_y
        ));
        const double y_out = 4148.0 * 0.075 - beam_x;
        CHECK(!mosflm::spot_on_panel(
            eiger, -beam_y, y_out, 0.0, true, beam_x, beam_y
        ));
        CHECK(!mosflm::spot_on_panel(
            eiger, -beam_y, -beam_x, 24.0, true, beam_x, beam_y
        ));
        const double x_in = 24.0 * 0.075 - beam_y;
        const double y_in = 24.0 * 0.075 - beam_x;
        CHECK(mosflm::spot_on_panel(
            eiger, x_in, y_in, 24.0, true, beam_x, beam_y
        ));
    }

    // Centered-beam fixtures keep the half-panel square.
    {
        DetectorGeometry sq;
        sq.nfast = 2000;
        sq.nslow = 2000;
        sq.pixel_size_x_mm = 0.1;
        sq.pixel_size_y_mm = 0.1;
        CHECK(mosflm::spot_on_panel(sq, 0.0, 0.0, 0.0, false, 0.0, 0.0));
        CHECK(mosflm::spot_on_panel(sq, 99.9, 0.0, 0.0, false, 0.0, 0.0));
        CHECK(!mosflm::spot_on_panel(sq, 100.0, 0.0, 0.0, false, 0.0, 0.0));
    }

    return mosflm_test::finish("detector");
}
