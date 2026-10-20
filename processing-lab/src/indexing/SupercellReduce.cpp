#include "mosflm/indexing/SupercellReduce.hpp"

#include <array>
#include <cmath>
#include <cstdlib>

namespace mosflm::indexing {
namespace {

/// Index-space directions m with entries in {-1,0,1} (one per ±m pair):
/// axes, face diagonals and body diagonals.
std::vector<Eigen::Vector3i> parity_directions() {
    std::vector<Eigen::Vector3i> out;
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            for (int k = -1; k <= 1; ++k) {
                const Eigen::Vector3i m(i, j, k);
                // Keep the first non-zero entry positive.
                const int lead = (i != 0) ? i : (j != 0) ? j : k;
                if (lead == 1) {
                    out.push_back(m);
                }
            }
        }
    }
    return out;
}

/// Integer basis of {h ∈ Z³ : m·h ≡ 0 (mod n)} as columns (det = ±n).
Eigen::Matrix3d sublattice_basis(const Eigen::Vector3i& m, int n) {
    int p = 0;
    while (m[p] == 0) {
        ++p;
    }
    Eigen::Matrix3d U = Eigen::Matrix3d::Zero();
    for (int q = 0; q < 3; ++q) {
        if (q == p) {
            U(p, q) = static_cast<double>(n);
            continue;
        }
        // e_q − m_p·m_q·e_p: m·(…) = m_q − m_p²·m_q = 0.
        U(q, q) = 1.0;
        U(p, q) = -static_cast<double>(m[p] * m[q]);
    }
    return U;
}

}  // namespace

SupercellReduceResult reduce_supercell_by_index_parity(
    const Eigen::Matrix3d& amat,
    const std::vector<float>& x,
    const std::vector<float>& y,
    const std::vector<float>& z,
    const SupercellReduceSettings& settings
) {
    SupercellReduceResult result;
    result.amat = amat;
    const std::vector<Eigen::Vector3i> dirs = parity_directions();
    constexpr std::array<int, 2> kFactors = {3, 2};

    for (int pass = 0; pass < 6; ++pass) {
        Eigen::FullPivLU<Eigen::Matrix3d> lu(result.amat);
        if (!lu.isInvertible()) {
            break;
        }
        const Eigen::Matrix3d Ainv = lu.inverse();
        std::vector<Eigen::Vector3i> hkls;
        hkls.reserve(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            const Eigen::Vector3d h =
                Ainv * Eigen::Vector3d(x[i], y[i], z[i]);
            const Eigen::Vector3d r = h.array().round();
            if ((h - r).cwiseAbs().maxCoeff() > settings.hkldev ||
                r.squaredNorm() < 0.5) {
                continue;
            }
            hkls.push_back(r.cast<int>());
        }
        if (static_cast<int>(hkls.size()) < settings.min_indexed) {
            break;
        }

        double best_frac = 0.0;
        int best_n = 0;
        Eigen::Vector3i best_m = Eigen::Vector3i::Zero();
        for (const int n : kFactors) {
            for (const auto& m : dirs) {
                int n_ok = 0;
                for (const auto& h : hkls) {
                    if (m.dot(h) % n == 0) {
                        ++n_ok;
                    }
                }
                const double frac = static_cast<double>(n_ok) /
                                    static_cast<double>(hkls.size());
                if (frac > best_frac + 1e-12) {
                    best_frac = frac;
                    best_n = n;
                    best_m = m;
                }
            }
        }
        if (best_n == 0 || best_frac < settings.min_frac) {
            break;
        }
        const Eigen::Matrix3d U = sublattice_basis(best_m, best_n);
        result.amat = result.amat * U;
        result.U = result.U * U;
        result.factor *= best_n;
    }
    return result;
}

}  // namespace mosflm::indexing
