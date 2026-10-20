/**
 * Kabsch Δ-vector basis recovery on a synthetic reciprocal lattice.
 */

#include "mosflm/indexing/DiffVectorIndex.hpp"
#include "test_util.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

using mosflm::indexing::DiffVectorSettings;
using mosflm::indexing::find_diff_vector_basis;

namespace {

bool near_pct(double got, double expect, double frac) {
    return std::fabs(got - expect) <= frac * std::fabs(expect);
}

}  // namespace

int main() {
    // Orthorhombic-ish direct cell ~60×80×100 Å → reciprocal axes.
    const double a = 60.0;
    const double b = 80.0;
    const double c = 100.0;
    const double wave = 1.0;
    Eigen::Matrix3d B = Eigen::Matrix3d::Zero();
    B(0, 0) = 1.0 / a;
    B(1, 1) = 1.0 / b;
    B(2, 2) = 1.0 / c;

    std::mt19937 rng(42);
    // Low-index reflections so primitive a*,b*,c* dominate the histogram.
    std::uniform_int_distribution<int> hkl(-5, 5);
    std::normal_distribution<double> noise(0.0, 0.00008);

    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> z;
    x.reserve(300);
    y.reserve(300);
    z.reserve(300);
    for (int n = 0; n < 300; ++n) {
        Eigen::Vector3d h(
            static_cast<double>(hkl(rng)),
            static_cast<double>(hkl(rng)),
            static_cast<double>(hkl(rng))
        );
        if (h.squaredNorm() < 1e-12) {
            continue;
        }
        Eigen::Vector3d r = B * h;
        r.x() += noise(rng);
        r.y() += noise(rng);
        r.z() += noise(rng);
        x.push_back(static_cast<float>(r.x()));
        y.push_back(static_cast<float>(r.y()));
        z.push_back(static_cast<float>(r.z()));
    }
    CHECK(x.size() >= 100u);

    DiffVectorSettings settings;
    settings.sepmin = 6.0f;
    settings.cluster_radius = 2.0f;
    settings.max_spots = 300;
    settings.max_diff_vectors = 100000;
    const auto result =
        find_diff_vector_basis(x, y, z, wave, settings);

    std::cout << "diff_vector_index: ok=" << result.ok
              << " diffs=" << result.n_diff_vectors
              << " clusters=" << result.n_clusters
              << " explained=" << result.n_explained << " cell=" << result.cell.a
              << " " << result.cell.b << " " << result.cell.c << '\n';
    CHECK(result.ok);
    CHECK(result.n_clusters >= 3u);
    CHECK(result.n_explained >= 20);

    // Edges may permute; match against {60,80,100}.
    std::vector<double> edges = {result.cell.a, result.cell.b, result.cell.c};
    std::sort(edges.begin(), edges.end());
    CHECK_MSG(near_pct(edges[0], 60.0, 0.12), "short edge not ~60");
    CHECK_MSG(near_pct(edges[1], 80.0, 0.12), "mid edge not ~80");
    CHECK_MSG(near_pct(edges[2], 100.0, 0.12), "long edge not ~100");

    return mosflm_test::finish("diff_vector_index");
}
