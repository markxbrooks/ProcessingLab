#pragma once

#include "mosflm/diffraction/MeasurementBox.hpp"
#include "mosflm/diffraction/SummationIntegrator.hpp"

#include <memory>
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
    /// MOSFLM PKRATIO: ``√(Σ (od − model)² / max(1, |bg| + |I_sum|))`` over
    /// valid peak pixels. About 1 when the profile matches counting
    /// statistics. ``FLAG`` bit 2 is set above 5 (MOSFLM's default PKRAT).
    double profile_gof{0.0};
    bool profile_ok{false};
    /// Detector-region index ``iy * n_bins_x + ix`` (0 if global-only).
    int profile_region{0};
};

/// Learned standard-profile sums: unit-volume learner profiles, sampled by
/// each pixel centre's offset from the predicted spot centre on an
/// oversampled grid, summed per detector region and over all regions, with
/// counts. Additive, so several frames' learners can be pooled into the
/// profiles for the next block.
struct ProfileAccumulator {
    int half{-1};  ///< Peak half-width the profiles were learned at.
    int oversample{1};  ///< Profile cells per pixel (each axis).
    int n_bins_x{0};
    int n_bins_y{0};
    /// [slot][cell] bilinear-splatted learner values and their weights; the
    /// profile is ``sum / w`` per cell. Slots are the ``n_bins_x × n_bins_y`` tiles,
    /// then six 8°-wide rings of 2θ (a tile's fallback).
    std::vector<double> region_sum;
    std::vector<double> region_w;
    std::vector<int> region_n;
    std::vector<double> global_sum;
    std::vector<double> global_w;
    int global_n{0};

    /// Add ``other``; ignored if its peak size or tiling differs.
    void add(const ProfileAccumulator& other);
    [[nodiscard]] bool empty() const { return global_n == 0; }
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
    /// If true, write profile intensity as the primary ``intensity`` field
    /// when ``pkratio`` is below ``max_pkratio_for_profile``.
    bool prefer_profile{true};
    /// Prefer summation when ``|I_sum−I_prf|/|I_prf|`` exceeds this.
    double max_pkratio_for_profile{1.5};
    /// Rectangular detector tiling (MOSFLM ``NMASKS`` / XLINES×YLINES).
    /// Defaults ``1×1`` (single global profile). PROCESS uses ``5×5`` (25).
    int n_bins_x{1};
    int n_bins_y{1};
    /// Learners pooled from earlier frames, added to this frame's before
    /// profiles are formed (MOSFLM / XDS average profiles over a block of
    /// images; one 0.15° Eiger frame left 15 of 25 regions short of
    /// learners, so 60% of spots were fitted with the global profile).
    std::shared_ptr<const ProfileAccumulator> prior{};
    /// Profile cells per pixel (each axis). Spots are learned and fitted at
    /// their predicted sub-pixel positions (XDS reflection frame); 1 is the
    /// old integer-pixel profile.
    int oversample{4};
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

    /// As above; also returns this frame's own learners (without the
    /// prior) in ``learned`` so callers can pool them.
    [[nodiscard]] std::vector<ProfileObservation> fit(
        const IntegrationRequest& request,
        const StandardProfileSettings& settings,
        ProfileAccumulator* learned
    ) const;
};

}  // namespace mosflm::diffraction
