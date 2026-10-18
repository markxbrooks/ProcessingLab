#pragma once

#include "mosflm/diffraction/MeasurementBox.hpp"
#include "mosflm/diffraction/SummationIntegrator.hpp"

#include <vector>

namespace mosflm::diffraction {

/// One standard-profile fitted observation (MOSFLM INTEG3/INTEG4 contract).
///
/// Intensity uses the unweighted Rossmann formula from ``unoptimizable.f``
/// PROCESS::
/// ``I = (Σ P·(OD − bg) · ΣP) / ΣP²`` with planar ``bg = a·p + b·q + c``.
struct ProfileObservation : MeasObservation {
    double profile_intensity{};
    double profile_sigma{};
    double summation_intensity{};
    double summation_sigma{};
    double pkratio{};  ///< ``|I_sum − I_prf| / max(|I_prf|,1)``
    bool profile_ok{false};
    /// Detector-region index ``iy * n_bins_x + ix`` (0 if global-only).
    int profile_region{0};
};

/// Options for building / applying detector-region standard profiles.
struct StandardProfileSettings {
    MeasurementBoxSettings box{};
    /// Minimum summation I/σ to contribute to the learned profile.
    double learn_min_isigi{5.0};
    /// Soft floor on number of learners per region; if fewer, fall back to
    /// the global (all-region) profile, then a Gaussian seed.
    int min_learners{8};
    /// Seed Gaussian σ (px) when the learned profile is too weak.
    double fallback_sigma_px{1.5};
    /// If true, write profile intensity as the primary ``intensity`` field.
    bool prefer_profile{true};
    /// Rectangular detector tiling (MOSFLM ``NMASKS`` / XLINES×YLINES).
    /// Defaults ``1×1`` (single global profile). PROCESS uses ``5×5`` (25).
    int n_bins_x{1};
    int n_bins_y{1};
};

/// Learn peak-shaped standard profiles per detector region, then fit spots.
///
/// Replaces the 2D-Gaussian :class:`ProfileFitter` stand-in for the
/// known-indexed PROCESS path. Regions are a uniform pixel grid
/// (``n_bins_x × n_bins_y``); spots with too few local learners use the
/// global accumulated profile.
class StandardProfileFitter {
public:
    [[nodiscard]] std::vector<ProfileObservation> fit(
        const IntegrationRequest& request,
        const StandardProfileSettings& settings
    ) const;
};

}  // namespace mosflm::diffraction
