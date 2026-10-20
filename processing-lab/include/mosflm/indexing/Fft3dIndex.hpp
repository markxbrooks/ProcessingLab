#pragma once

#include "mosflm/geometry/UnitCell.hpp"

#include <Eigen/Dense>

#include <cstddef>
#include <vector>

namespace mosflm::indexing {

/// DIALS-style 3D-FFT lattice search (native C++; no CCTBX/DIALS link).
///
/// Bins reciprocal lattice points onto a cubic grid, takes a 3D FFT, and
/// interprets magnitude peaks as candidate direct-space basis vectors.
struct Fft3dSettings {
    /// Grid edge length (power of two, clamped to ``[32, 128]``).
    /// Auto-bumped to 128 when ``dmax`` / ``b_max`` exceeds 100 Å.
    int grid_n{128};
    /// Cap input RLVs (histogram fill).
    std::size_t max_spots{8000};
    /// Keep this many FFT peaks before basis search.
    std::size_t n_peaks{40};
    /// Minimum |direct| peak length (Å); suppress origin neighbours.
    float b_min_A{15.0f};
    /// Maximum |direct| peak length (Å); 0 → use ``dmax_A``.
    float b_max_A{0.0f};
    /// Max |hkl − round(hkl)| when scoring a candidate basis.
    float index_tol{0.30f};
};

struct Fft3dResult {
    bool ok{false};
    /// Reciprocal orientation matrix for ``BravaisClassifier`` (row i = basis i).
    Eigen::Matrix3d reciprocal_om{Eigen::Matrix3d::Identity()};
    UnitCell cell{};
    int grid_n{0};
    std::size_t n_spots_used{0};
    std::size_t n_peaks{0};
    int n_indexed{0};
};

/// Find a reciprocal lattice basis via 3D FFT of reciprocal ``xyz`` (Å⁻¹).
///
/// \param x,y,z Reciprocal coordinates in the same frame as AutoIndexer DPS.
/// \param wavelength_A Unused for sizing (kept for API symmetry with Δ-vector).
/// \param dmax_A Largest expected direct-cell edge (Å); sets real-space search.
[[nodiscard]] Fft3dResult find_fft3d_basis(
    const std::vector<float>& x,
    const std::vector<float>& y,
    const std::vector<float>& z,
    double wavelength_A,
    const Fft3dSettings& settings = {},
    double dmax_A = 200.0
);

}  // namespace mosflm::indexing
