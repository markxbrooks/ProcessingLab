#pragma once

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/crystal/MtzAdapter.hpp"
#include "mosflm/diffraction/MeasurementBox.hpp"
#include "mosflm/diffraction/Shoebox.hpp"
#include "mosflm/diffraction/SummationIntegrator.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"

#include <vector>

namespace mosflm::diffraction {

/// Controls for 3D standard-profile integration (XDS/DIALS-style shoebox).
struct Profile3dSettings {
    MeasurementBoxSettings box{};
    /// Shoebox depth in frames (odd; even values round up). CONTROL ``NADD``.
    int nadd{5};
    double learn_min_isigi{3.0};
    int min_learners{8};
    double fallback_sigma_xy_px{1.5};
    double fallback_sigma_phi_frames{1.0};
    int n_bins_x{5};
    int n_bins_y{5};
    /// Extra constant-φ term on the background plane (``bg += d·Δφ``).
    bool use_phi_bg{true};
    bool prefer_profile{true};
    double max_pkratio_for_profile{2.0};
    /// Mosaic / oscillation for fallback φ Gaussian and FRACTIONCALC.
    double mosaic_spread_deg{0.15};
    double osc_deg{0.15};
};

/// Planar + optional φ background: ``bg = a·Δcol + b·Δrow + c + d·Δφ``.
struct Background3d {
    double a{0.0};
    double b{0.0};
    double c{0.0};
    double d{0.0};
    bool ok{false};
    int n_used{0};

    [[nodiscard]] double at(double dc, double dr, double dphi) const {
        return a * dc + b * dr + c + d * dphi;
    }
};

/// One 3D-profile observation (one row per HKL).
struct Profile3dObservation {
    int h{};
    int k{};
    int l{};
    double x_mm{};
    double y_mm{};
    double phi_calc_deg{};
    int peak_frame{};
    int col0{};
    int row0{};
    double intensity{};
    double sigma{};
    double summation_intensity{};
    double summation_sigma{};
    double profile_intensity{};
    double profile_sigma{};
    double fractioncalc{1.0};
    double pkratio{0.0};
    int npix_peak{};
    int npix_bg{};
    int nphi{};
    int profile_region{-1};
    bool ok{false};
    bool profile_ok{false};
};

/// Learn a 3D Rossmann profile per NMASKS tile and fit each shoebox once.
class Profile3dFitter {
public:
    [[nodiscard]] std::vector<Profile3dObservation> fit(
        const std::vector<Shoebox>& boxes,
        const DetectorGeometry& detector,
        const IntegrationSettings& integration,
        const Profile3dSettings& settings
    ) const;
};

/// Convert a 3D observation to a MOSFLM PROCESS MTZ row.
///
/// Does **not** scale by ``1/FRACTIONCALC``: the shoebox already covers the
/// rocking curve. Applies rotation-method LP when ``integration.rotation_lp``.
[[nodiscard]] crystal::MosflmReflection profile3d_to_mosflm(
    const Profile3dObservation& obs,
    const DetectorGeometry& detector,
    const OrientationMatrix& orientation,
    const crystal::CrystalMetadata& crystal,
    const IntegrationSettings& integration,
    double air_absorption_mm_inv = 0.0,
    double osc_deg = 0.0
);

}  // namespace mosflm::diffraction
