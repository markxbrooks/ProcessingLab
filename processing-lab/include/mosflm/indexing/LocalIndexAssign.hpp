#pragma once

#include <Eigen/Dense>

#include <cstddef>
#include <vector>

namespace mosflm::indexing {

/// Settings for Kabsch-style local indexing (XDS IDXREF §2.7).
struct LocalIndexSettings {
    /// Max |fractional Δhkl − round| for a reliable branch (Kabsch ε).
    double epsilon{0.05};
    /// Max |integer Δhkl| on a branch (Kabsch δ).
    int delta_max{5};
    /// Branches longer than this start a new subtree (aliens).
    double ell_min{0.05};
    /// Cap spots for MST (O(n²) edges); brightest / first n used.
    std::size_t max_spots{800};
};

/// Per-spot Miller assignment from local indexing.
struct LocalIndexSpot {
    int h{0};
    int k{0};
    int l{0};
    bool indexed{false};
    int subtree{0};
};

struct LocalIndexResult {
    std::vector<LocalIndexSpot> spots;
    std::size_t n_indexed{0};
    std::size_t n_largest_subtree{0};
};

/**
 * Assign HKL by a shortest spanning tree on reciprocal lattice points.
 *
 * Edge length follows Kabsch: short when Δr ≈ A·Δhkl with small integers.
 * Walks from the densest hub; spots in small alien subtrees are unmarked.
 *
 * \param amat MOSFLM orientation (x = A · hkl).
 * \param x,y,z Reciprocal lab coordinates (Å⁻¹), same frame as AutoIndexer.
 */
[[nodiscard]] LocalIndexResult assign_local_indices(
    const Eigen::Matrix3d& amat,
    const std::vector<float>& x,
    const std::vector<float>& y,
    const std::vector<float>& z,
    const LocalIndexSettings& settings = {}
);

}  // namespace mosflm::indexing
