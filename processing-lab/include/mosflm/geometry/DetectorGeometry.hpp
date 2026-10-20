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
    /// MOSFLM YSCALE: slow-axis pixel aspect (``mm_to_pixel``); 1 = square.
    double yscale{1.0};
    /// MOSFLM refined distance scale (``XTOFRA``); multiplies ideal mm.
    double xtofra{1.0};
    /// Detector omega (``CCOM``) in degrees; rotates ideal XG/YG in MMTOPX.
    double ccom_deg{0.0};
    /// Detector tilt / twist as printed ITILT/ITWIST (MOSFLM milliradian ints).
    double tilt_mrad{0.0};
    double twist_mrad{0.0};
    /// When true, use ``ROTMAT2``/``NEWT`` ray–plane scale; else old MMTOPX
    /// ``DXTOF = XTOFRA + TILT·X + TWIST·Y`` (default, matches insulin oracle).
    bool use_rotmat2{false};

    /**
     * Detector mm for a scattering vector that intersects the detector.
     *
     * Flat projection at ``z = distance`` (beam along +z), then MOSFLM
     * ``MMTOPX`` tilt/twist distance scale in the XG/YG frame
     * (``XG↔y_mm``, ``YG↔x_mm``). YSCALE is applied in :meth:`mm_to_pixel`.
     *
     * Returns false if the ray does not hit the plane with positive z.
     */
    [[nodiscard]] bool lab_to_detector_mm(
        const Eigen::Vector3d& s1,
        double& x_mm,
        double& y_mm
    ) const;

    void mm_to_pixel(double x_mm, double y_mm, double& xf, double& yf) const;

    [[nodiscard]] double two_theta_deg(double x_mm, double y_mm) const;
};

/**
 * Highest resolution (smallest d, Å) among the four panel corners.
 *
 * Corners are inset by ``margin_px``. ``beam_x_mm`` / ``beam_y_mm`` are the
 * beam in millimetres from the panel origin. PROCESS leaves
 * ``DetectorGeometry::beam_*`` at 0, so the beam is passed separately.
 * The result is reduced by a relative 1e-4 so ``d < dmin`` keeps a reflection
 * that lands on that pixel. Returns 0 when the geometry cannot be used.
 */
[[nodiscard]] double detector_corner_dmin_A(
    const DetectorGeometry& det,
    double beam_x_mm,
    double beam_y_mm,
    double margin_px
);

/**
 * True when a beam-relative C++ detector position (mm) lies on the panel.
 *
 * ``use_panel_beam`` false: square centered on the beam,
 * ``|pixel| < n/2 - margin`` (legacy fixtures, beam at the panel centre).
 * ``use_panel_beam`` true: MOSFLM ``XDET``/``YDET`` inside
 * ``[margin, n - margin)``.
 */
[[nodiscard]] bool spot_on_panel(
    const DetectorGeometry& det,
    double x_mm,
    double y_mm,
    double margin_px,
    bool use_panel_beam,
    double panel_beam_x_mm,
    double panel_beam_y_mm
);

/**
 * MOSFLM ``ROTMAT2``: detector normal from TILT/TWIST (radians).
 *
 * ``DETNOR = (sin(-TILT), sin(-TWIST)·cos(-TILT), cos(-TWIST)·cos(-TILT))``
 * with flat detector ``(0,0,1)``.
 */
[[nodiscard]] Eigen::Vector3d rotmat2_detnor(double tilt_rad, double twist_rad);

}  // namespace mosflm
