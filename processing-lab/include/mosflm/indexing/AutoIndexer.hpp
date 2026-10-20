#pragma once

#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/geometry/UnitCell.hpp"
#include "mosflm/indexing/BravaisClassifier.hpp"
#include "mosflm/indexing/SpotList.hpp"

#include <Eigen/Dense>

#include <vector>

namespace mosflm::indexing {

/// Detector / beam / wavelength for converting ``.spt`` spots to reciprocal xyz.
struct AutoIndexGeometry {
    double distance_mm{160.0};
    double wavelength_A{0.9790};
    double beam_x_mm{94.0};
    double beam_y_mm{94.7};
    double pixel_mm{0.0816};
    double yscale{1.0};
    /// Maximum expected cell edge (Å); MOSFLM ``DMAX``.
    double dmax_A{200.0};
};

/// Knobs matching MOSFLM ``/DPSDATA/`` defaults (+ insulin ``ipmosflm.in``).
struct AutoIndexSettings {
    int nvect{30};
    int nvref{30};
    int maxcyc{10};
    int inicyc{3};
    float convrg{0.0004f};
    int numormat{50};
    float stepscal{1.0f};
    int numrlp{2000};
    int numsptdps{2000};
    int npkfftmin{6};
    int nmultfft{2};
    float theta1d{0.0f};
    float theta2d{0.0f};
    float sinthmin{0.5f};
    float volmin{0.5f};
    float hkldev{0.3f};
    float hklmin{0.75f};
    float dmaxscal{1.5f};
    /// Keep spots with ``I/σ ≥ min_isigi`` (MOSFLM AUTOINDEX THRESHOLD).
    float min_isigi{20.0f};
    /// Prefer lattice character (1–44). 0 = SUGGLATT / prefer_system.
    /// Insulin solution 12 maps to a cI-compatible character after Bravais.
    int prefer_character{0};
    /// If non-empty (e.g. ``"cI"``), pick lowest-criterion matching system.
    std::string prefer_system;
    /// Max Bravais penalty for SUGGLATT-style highest-symmetry pick
    /// (DiffractionLab ``max_penalty_to_refine`` / Mosflm soft cut).
    float max_penalty{50.0f};
    /// BESTSOLAR-lite / XDS IDXREF: require at least this fraction of spots
    /// to index (``|Δh|,|Δk|,|Δl| ≤ hkldev``). Default 0.70 matches XDS.
    float min_indexed_frac{0.70f};
    /// Max RMS fractional HKL residual for indexed spots.
    float max_rmsd_hkl{0.25f};
    /// Kabsch Δ-vector lattice search (XDS IDXREF) as a second candidate.
    bool use_diff_vector{true};
    /// DIALS-style native 3D-FFT lattice search (third INDEXFIT candidate).
    bool use_fft3d{true};
    /// Run Cambridge DPS ``index_fortran_`` (disable with CONTROL ``FFTONLY``).
    bool use_dps{true};
    /// Geometry LS after INDEXFIT (Kabsch §2.8-lite). Off with ``NOIDXREF``.
    bool use_idxref_refine{true};
    /// XDS ``SEPMIN=`` (units of ``0.001/WAVE``).
    float sepmin{6.0f};
    /// XDS ``CLUSTER_RADIUS=`` (same units as SEPMIN).
    float cluster_radius{2.0f};
};

/// Spot-vs-AMAT indexing quality (BESTSOLAR-lite).
struct IndexFitStats {
    double indexed_frac{0.0};
    double rmsd_hkl{0.0};
    std::size_t n_indexed{0};
    std::size_t n_total{0};
};

struct AutoIndexResult {
    bool ok{false};
    int ierr{0};
    UnitCell cell{};
    /// MOSFLM Orientation Matrix [A] (rows = lab X,Y,Z).
    Eigen::Matrix3d amat{Eigen::Matrix3d::Identity()};
    OrientationMatrix orientation{};
    std::vector<BravaisCandidate> bravais;
    int selected_character{0};
    std::string selected_system;
    /// Default space-group number for ``selected_system`` (Mosflm SPGMIN).
    int selected_symmetry_number{1};
    std::size_t n_spots_used{0};
    IndexFitStats fit{};
    /// ``"dps"``, ``"delta"``, or ``"fft3d"`` — INDEXFIT winner.
    std::string index_method{"dps"};
    /// IDXREF geometry LS was applied: ``amat``/``cell`` are only consistent
    /// with these beam / distance values (``.spt`` mm, same frame as input).
    bool idxref_applied{false};
    double idxref_beam_x_mm{0.0};
    double idxref_beam_y_mm{0.0};
    double idxref_distance_mm{0.0};
};

/// Stage 1 AUTOINDEX façade: MOSFLM ``.spt`` → DPS ``index_fortran_`` → AMAT.
///
/// Links existing DPS archives (``libdps_index.a`` / ``libdps.a`` /
/// ``dps_index.a``). Does not reimplement FFT indexing or spot finding.
class AutoIndexer {
public:
    [[nodiscard]] AutoIndexResult index(
        const SpotList& spots,
        const AutoIndexGeometry& geometry,
        const AutoIndexSettings& settings = {}
    ) const;
};

}  // namespace mosflm::indexing
