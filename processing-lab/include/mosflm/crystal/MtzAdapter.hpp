#pragma once

#include "mosflm/crystal/CrystalMetadata.hpp"

#include <limits>
#include <filesystem>
#include <string>
#include <vector>

namespace mosflm::crystal {

struct SimpleReflection {
    int h{};
    int k{};
    int l{};
    float intensity{};
    float sigma{};
};

/// One MOSFLM-style PROCESS MTZ row (STARTMTZ / WRMTZ ``MCOLS=18``).
struct MosflmReflection {
    int h{};
    int k{};
    int l{};
    int misym{1};  ///< packed ``M/ISYM`` (``MP*256 + KSYM``)
    int batch{1};
    float intensity{};
    float sigma{};
    float ipr{};      ///< profile intensity
    float sigipr{};
    float fractioncalc{1.0f};
    float xdet{};
    float ydet{};
    float rot{};
    float width{};
    float lp{1.0f};
    int mpart{0};
    int flag{0};  ///< 0 = good; non-zero = IBAD-style bits (Aimless rejects)
    float bgpkratios{0.0f};
    /// Observed centroid (XDET/YDET pixels) for geometry refinement; NaN if
    /// not measured. In-memory only — not written to the MTZ.
    float xdet_obs{std::numeric_limits<float>::quiet_NaN()};
    float ydet_obs{std::numeric_limits<float>::quiet_NaN()};
    /// Original (pre-ASU) Miller indices, for partial summation and
    /// refinement. In-memory only.
    int h0{};
    int k0{};
    int l0{};
    /// Observed φ (degrees): frame centre for a slice; intensity-weighted
    /// slice φ after :func:`sum_partial_slices`. NaN if unknown.
    float phi_obs{std::numeric_limits<float>::quiet_NaN()};
    /// Intensity-weighted SD of slice φ about ``phi_obs`` (degrees), set by
    /// :func:`sum_partial_slices` when ≥2 slices carry signal.
    float phi_obs_sd{std::numeric_limits<float>::quiet_NaN()};
};

/// Rotation range of one image, for MTZ batch headers.
struct BatchScan {
    int batch{};
    double phi_start_deg{};
    double phi_end_deg{};
};

/// Minimal MTZ writer/reader for integrate comparisons (Gemmi Mtz API).
class MtzAdapter {
public:
    /// Write a single-dataset MTZ with H K L I SIGI columns.
    static void write_simple(
        const std::filesystem::path& path,
        const CrystalMetadata& crystal,
        const std::vector<SimpleReflection>& reflections
    );

    /// Write a MOSFLM-compatible 18-column PROCESS MTZ.
    ///
    /// Emits one MTZ batch header per unique ``BATCH`` value so CCP4 tools
    /// (POINTLESS / AIMLESS) treat the file as unmerged. ``wavelength_A`` is
    /// written to the dataset and each batch header when ``> 0``. ``scans``
    /// fills each batch's φ start / end / range (AIMLESS uses them).
    static void write_mosflm(
        const std::filesystem::path& path,
        const CrystalMetadata& crystal,
        const std::vector<MosflmReflection>& reflections,
        double wavelength_A = 0.0,
        const std::vector<BatchScan>& scans = {},
        bool reindex_to_asu = false
    );

    /// Read H K L I SIGI from an MTZ (works for simple or full MOSFLM files).
    [[nodiscard]] static std::vector<SimpleReflection> read_simple(
        const std::filesystem::path& path
    );

    /// Reduce ``(h,k,l)`` to the space-group ASU and return packed ``M/ISYM``.
    ///
    /// Uses Gemmi ``ReciprocalAsu``; ``misym = isym`` (compatible with AIMLESS
    /// when ``M=0`` for centrics / standard packing).
    static void reduce_to_asu(
        const gemmi::SpaceGroup* sg,
        int& h,
        int& k,
        int& l,
        int& misym
    );
};

}  // namespace mosflm::crystal
