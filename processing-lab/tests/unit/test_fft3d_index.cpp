/**
 * Native FFT3D basis recovery: orthorhombic + tetragonal (√2 supercell trap).
 */

#include "mosflm/indexing/Fft3dIndex.hpp"
#include "test_util.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

using mosflm::indexing::Fft3dSettings;
using mosflm::indexing::find_fft3d_basis;

namespace {

bool near_pct(double got, double expect, double frac) {
    return std::fabs(got - expect) <= frac * std::fabs(expect);
}

void fill_rlps(
    const Eigen::Matrix3d& B,
    int n_try,
    unsigned seed,
    std::vector<float>& x,
    std::vector<float>& y,
    std::vector<float>& z
) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> hkl(-6, 6);
    std::normal_distribution<double> noise(0.0, 0.00005);
    x.clear();
    y.clear();
    z.clear();
    x.reserve(static_cast<std::size_t>(n_try));
    y.reserve(static_cast<std::size_t>(n_try));
    z.reserve(static_cast<std::size_t>(n_try));
    for (int n = 0; n < n_try; ++n) {
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
}

}  // namespace

int main() {
    const double wave = 1.0;

    // --- Orthorhombic 60×80×100 ---
    {
        Eigen::Matrix3d B = Eigen::Matrix3d::Zero();
        B(0, 0) = 1.0 / 60.0;
        B(1, 1) = 1.0 / 80.0;
        B(2, 2) = 1.0 / 100.0;
        std::vector<float> x;
        std::vector<float> y;
        std::vector<float> z;
        fill_rlps(B, 500, 7, x, y, z);
        CHECK(x.size() >= 200u);

        Fft3dSettings settings;
        settings.grid_n = 128;
        settings.max_spots = 500;
        settings.n_peaks = 48;
        settings.b_min_A = 20.0f;
        settings.b_max_A = 120.0f;
        const auto result =
            find_fft3d_basis(x, y, z, wave, settings, /*dmax_A=*/120.0);

        std::cout << "fft3d_index ortho: ok=" << result.ok
                  << " cell=" << result.cell.a << " " << result.cell.b << " "
                  << result.cell.c << '\n';
        CHECK(result.ok);
        CHECK(result.n_indexed >= 40);
        std::vector<double> edges = {
            result.cell.a, result.cell.b, result.cell.c
        };
        std::sort(edges.begin(), edges.end());
        CHECK_MSG(near_pct(edges[0], 60.0, 0.08), "ortho short not ~60");
        CHECK_MSG(near_pct(edges[1], 80.0, 0.08), "ortho mid not ~80");
        CHECK_MSG(near_pct(edges[2], 100.0, 0.08), "ortho long not ~100");
    }

    // --- Tetragonal 58×58×150 (thaumatin-like; reject ~82×82 supercell) ---
    {
        Eigen::Matrix3d B = Eigen::Matrix3d::Zero();
        B(0, 0) = 1.0 / 58.0;
        B(1, 1) = 1.0 / 58.0;
        B(2, 2) = 1.0 / 150.0;
        std::vector<float> x;
        std::vector<float> y;
        std::vector<float> z;
        fill_rlps(B, 800, 42, x, y, z);
        CHECK(x.size() >= 300u);

        Fft3dSettings settings;
        settings.grid_n = 128;
        settings.max_spots = 800;
        settings.n_peaks = 48;
        settings.b_min_A = 20.0f;
        settings.b_max_A = 200.0f;
        const auto result =
            find_fft3d_basis(x, y, z, wave, settings, /*dmax_A=*/200.0);

        std::cout << "fft3d_index tet: ok=" << result.ok
                  << " cell=" << result.cell.a << " " << result.cell.b << " "
                  << result.cell.c << " indexed=" << result.n_indexed << '\n';
        CHECK(result.ok);
        CHECK(result.n_indexed >= 40);

        std::vector<double> edges = {
            result.cell.a, result.cell.b, result.cell.c
        };
        std::sort(edges.begin(), edges.end());
        // Basal edges ~58 (not ~58√2≈82); long axis ~150.
        CHECK_MSG(near_pct(edges[0], 58.0, 0.08), "tet short basal not ~58");
        CHECK_MSG(near_pct(edges[1], 58.0, 0.08), "tet other basal not ~58");
        CHECK_MSG(near_pct(edges[2], 150.0, 0.08), "tet c not ~150");
        CHECK_MSG(
            !near_pct(edges[0], 58.0 * std::sqrt(2.0), 0.06),
            "tet locked onto √2 supercell ~82"
        );
    }

    return mosflm_test::finish("fft3d_index");
}
