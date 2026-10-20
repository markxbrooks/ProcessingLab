#pragma once

#include "mosflm/geometry/UnitCell.hpp"

#include <Eigen/Dense>

#include <cstddef>
#include <vector>

namespace mosflm::indexing {

/// XDS IDXREF-style reciprocal difference-vector lattice search (Kabsch).
///
/// ``SEPMIN`` / ``CLUSTER_RADIUS`` use XDS units of ``0.001 / wavelength``
/// (Å⁻¹). Defaults match XDS.INP (``SEPMIN=6``, ``CLUSTER_RADIUS=2``).
struct DiffVectorSettings {
    float sepmin{6.0f};
    float cluster_radius{2.0f};
    /// Cap input RLVs (pair enumeration is O(n²)).
    std::size_t max_spots{500};
    /// Hard cap on accepted difference vectors (keeps clustering O(n)).
    std::size_t max_diff_vectors{40000};
    /// Keep the densest clusters for basis search.
    std::size_t max_clusters{60};
    /// Try triplets among the top-N densest clusters.
    std::size_t basis_pool{18};
    /// Max |hkl − round(hkl)| to count a cluster as explained by a basis.
    float index_tol{0.25f};
};

struct DiffVectorCluster {
    Eigen::Vector3d center{Eigen::Vector3d::Zero()};
    int count{0};
};

struct DiffVectorResult {
    bool ok{false};
    /// Reciprocal orientation matrix for ``BravaisClassifier`` (row i = basis i).
    Eigen::Matrix3d reciprocal_om{Eigen::Matrix3d::Identity()};
    UnitCell cell{};
    std::size_t n_diff_vectors{0};
    std::size_t n_clusters{0};
    int n_explained{0};
    std::vector<DiffVectorCluster> clusters;
};

/// Find a reciprocal lattice basis from Mosflm/DPS reciprocal ``xyz`` (Å⁻¹).
///
/// \param x,y,z Reciprocal coordinates in the same frame as AutoIndexer DPS.
/// \param wavelength_A Beam wavelength (Å); sets the XDS length unit.
/// \param dmax_A Largest expected direct-cell edge (Å). When ``> 0``, SEPMIN
///        is capped so that ``|b*| ≈ 1/dmax`` is not excluded (large cells).
[[nodiscard]] DiffVectorResult find_diff_vector_basis(
    const std::vector<float>& x,
    const std::vector<float>& y,
    const std::vector<float>& z,
    double wavelength_A,
    const DiffVectorSettings& settings = {},
    double dmax_A = 0.0
);

}  // namespace mosflm::indexing
