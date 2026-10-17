#pragma once

#include <Eigen/Dense>
#include <cmath>

namespace mosflm {

/// Flat detector geometry in millimetres / pixels (MOSFLM-style).
struct DetectorGeometry {
    double distance_mm{100.0};
    double beam_x_mm{0.0};  ///< Beam centre relative to detector origin
    double beam_y_mm{0.0};
    double pixel_size_x_mm{0.1};
    double pixel_size_y_mm{0.1};
    int nfast{1000};
    int nslow{1000};
    double wavelength_A{1.0};

    /// Detector mm coordinates for a scattering vector that intersects the
    /// detector plane at z = distance (simple flat detector, beam along +z).
    /// Returns false if the ray does not hit the plane with positive z.
    [[nodiscard]] bool lab_to_detector_mm(
        const Eigen::Vector3d& s1,
        double& x_mm,
        double& y_mm
    ) const;

    void mm_to_pixel(double x_mm, double y_mm, double& xf, double& yf) const;

    [[nodiscard]] double two_theta_deg(double x_mm, double y_mm) const;
};

}  // namespace mosflm
