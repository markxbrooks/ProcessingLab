#include "mosflm/indexing/DiffVectorIndex.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace mosflm::indexing {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kRad2Deg = 180.0 / kPi;

/// Direct cell from reciprocal basis rows of ``om`` (Å⁻¹).
UnitCell cell_from_reciprocal_basis(const Eigen::Matrix3d& om) {
    const Eigen::Vector3d b0 = om.row(0);
    const Eigen::Vector3d b1 = om.row(1);
    const Eigen::Vector3d b2 = om.row(2);
    const double vst = b0.dot(b1.cross(b2));
    UnitCell out;
    if (!(std::fabs(vst) > 1e-18)) {
        return out;
    }
    const double inv_v = 1.0 / std::fabs(vst);
    out.a = b1.cross(b2).norm() * inv_v;
    out.b = b2.cross(b0).norm() * inv_v;
    out.c = b0.cross(b1).norm() * inv_v;
    // Reciprocal metric → direct angles via cos α = (b*·c*)/(|b||c|) wait:
    // Direct α is angle between b and c real-space axes = angle between
    // (b2×b0) and (b0×b1) normals... Use G = (B^T B)^{-1} with B cols = bi*.
    Eigen::Matrix3d B;
    B.col(0) = b0;
    B.col(1) = b1;
    B.col(2) = b2;
    Eigen::Matrix3d Gstar = B.transpose() * B;
    Eigen::FullPivLU<Eigen::Matrix3d> lu(Gstar);
    if (!lu.isInvertible()) {
        return out;
    }
    const Eigen::Matrix3d G = lu.inverse();
    auto ang = [&](int i, int j, double li, double lj) {
        const double c = G(i, j) / (li * lj);
        const double cc = std::max(-1.0, std::min(1.0, c));
        return std::acos(cc) * kRad2Deg;
    };
    out.alpha_deg = ang(1, 2, out.b, out.c);
    out.beta_deg = ang(0, 2, out.a, out.c);
    out.gamma_deg = ang(0, 1, out.a, out.b);
    return out;
}

/// Flip so the largest-magnitude component is non-negative (hemisphere fold).
Eigen::Vector3d fold_hemisphere(Eigen::Vector3d v) {
    const double ax = std::fabs(v.x());
    const double ay = std::fabs(v.y());
    const double az = std::fabs(v.z());
    if (ax >= ay && ax >= az) {
        if (v.x() < 0.0) {
            v = -v;
        }
    } else if (ay >= az) {
        if (v.y() < 0.0) {
            v = -v;
        }
    } else if (v.z() < 0.0) {
        v = -v;
    }
    return v;
}

bool nearly_independent(
    const Eigen::Vector3d& a,
    const Eigen::Vector3d& b,
    const Eigen::Vector3d& c,
    double min_abs_det
) {
    const double det = a.dot(b.cross(c));
    return std::fabs(det) >= min_abs_det;
}

}  // namespace

DiffVectorResult find_diff_vector_basis(
    const std::vector<float>& x,
    const std::vector<float>& y,
    const std::vector<float>& z,
    double wavelength_A,
    const DiffVectorSettings& settings,
    double dmax_A
) {
    DiffVectorResult result;
    const std::size_t n_in = x.size();
    if (n_in < 20 || y.size() != n_in || z.size() != n_in) {
        return result;
    }
    if (!(wavelength_A > 1e-6)) {
        return result;
    }

    // XDS: length unit for SEPMIN / CLUSTER_RADIUS is 0.001 / WAVE (Å⁻¹).
    const double unit = 0.001 / wavelength_A;
    double sepmin = static_cast<double>(settings.sepmin) * unit;
    double crad = static_cast<double>(settings.cluster_radius) * unit;
    // Large cells: |b*| ≈ 1/dmax can fall below default SEPMIN=6 → exclude the
    // true basis. Cap physical SEPMIN at ~0.45/dmax (XDS large-cell advice).
    if (dmax_A > 50.0) {
        const double sep_cap = 0.45 / dmax_A;
        if (sep_cap < sepmin) {
            const double scale = sep_cap / sepmin;
            sepmin = sep_cap;
            crad *= scale;
        }
    }
    if (!(sepmin > 0.0) || !(crad > 0.0)) {
        return result;
    }
    // Reject very long differences (noise / unrelated spots).
    const double sepmax = std::max(sepmin * 40.0, 0.15);

    const std::size_t n =
        std::min(n_in, std::max<std::size_t>(20, settings.max_spots));

    std::vector<Eigen::Vector3d> rlv;
    rlv.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        rlv.emplace_back(
            static_cast<double>(x[i]),
            static_cast<double>(y[i]),
            static_cast<double>(z[i])
        );
    }

    std::vector<Eigen::Vector3d> diffs;
    diffs.reserve(n * 8);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            Eigen::Vector3d d = rlv[i] - rlv[j];
            const double len = d.norm();
            if (len < sepmin || len > sepmax) {
                continue;
            }
            diffs.push_back(fold_hemisphere(std::move(d)));
        }
    }
    result.n_diff_vectors = diffs.size();
    if (diffs.size() < 30) {
        return result;
    }

    // Greedy clusters within CLUSTER_RADIUS (stable on synthetic + real data).
    std::vector<DiffVectorCluster> clusters;
    clusters.reserve(512);
    for (const auto& d : diffs) {
        int best = -1;
        double best_d2 = crad * crad;
        for (std::size_t c = 0; c < clusters.size(); ++c) {
            const double d2 = (d - clusters[c].center).squaredNorm();
            if (d2 < best_d2) {
                best_d2 = d2;
                best = static_cast<int>(c);
            }
        }
        if (best < 0) {
            DiffVectorCluster cl;
            cl.center = d;
            cl.count = 1;
            clusters.push_back(cl);
        } else {
            auto& cl = clusters[static_cast<std::size_t>(best)];
            const double w = static_cast<double>(cl.count);
            cl.center = (cl.center * w + d) / (w + 1.0);
            ++cl.count;
        }
    }

    // Prefer dense then short clusters (short reciprocal ⇒ long direct edges).
    std::stable_sort(
        clusters.begin(), clusters.end(),
        [](const DiffVectorCluster& a, const DiffVectorCluster& b) {
            if (a.count != b.count) {
                return a.count > b.count;
            }
            return a.center.squaredNorm() < b.center.squaredNorm();
        }
    );
    // Keep a long list for basis shortening; rank/truncated a working set.
    std::vector<DiffVectorCluster> clusters_all = clusters;
    const std::size_t cluster_cap =
        std::max(settings.max_clusters, static_cast<std::size_t>(120));
    if (clusters.size() > cluster_cap) {
        clusters.resize(cluster_cap);
    }
    result.clusters = clusters;
    result.n_clusters = clusters.size();
    if (clusters.size() < 3) {
        return result;
    }

    // Expected RLV length window from DMAX (and a soft lower cell bound).
    const double len_lo =
        (dmax_A > 50.0) ? (0.35 / dmax_A) : (0.5 * sepmin);
    const double len_hi =
        (dmax_A > 50.0) ? std::max(3.0 / 25.0, 2.5 / 40.0)
                        : (sepmin * 25.0);

    // Basis pool: densest clusters in the expected |b*| window.
    const int min_count = std::max(
        3, static_cast<int>(std::lround(
               0.001 * static_cast<double>(result.n_diff_vectors)
           ))
    );
    const std::size_t pool_cap =
        std::max(settings.basis_pool, static_cast<std::size_t>(28));
    std::vector<std::size_t> pool_idx;
    pool_idx.reserve(pool_cap);
    for (std::size_t i = 0;
         i < clusters.size() && pool_idx.size() < pool_cap;
         ++i) {
        if (clusters[i].count < min_count) {
            continue;
        }
        const double len = clusters[i].center.norm();
        if (len < len_lo || len > len_hi) {
            continue;
        }
        pool_idx.push_back(i);
    }
    if (pool_idx.size() < 3) {
        pool_idx.clear();
        for (std::size_t i = 0;
             i < clusters.size() && pool_idx.size() < pool_cap;
             ++i) {
            if (clusters[i].center.norm() < len_lo) {
                continue;
            }
            pool_idx.push_back(i);
        }
    }
    if (pool_idx.size() < 3) {
        return result;
    }

    int total_weight = 0;
    for (const auto& cl : clusters) {
        total_weight += cl.count;
    }
    // Minimum |scalar triple product|; reject near-coplanar / tiny-volume bases.
    const double min_det = std::pow(sepmin, 3.0) * 0.5;

    double best_explained_frac = -1.0;
    double best_vol = -1.0;
    double best_basis_len = 1e9;
    double best_mean_idx = 1e9;
    Eigen::Matrix3d best_om = Eigen::Matrix3d::Identity();
    int best_explained = 0;

    auto score_basis = [&](const Eigen::Vector3d& b0,
                           const Eigen::Vector3d& b1,
                           const Eigen::Vector3d& b2) {
        Eigen::Matrix3d B;
        B.col(0) = b0;
        B.col(1) = b1;
        B.col(2) = b2;
        Eigen::FullPivLU<Eigen::Matrix3d> lu(B);
        if (!lu.isInvertible()) {
            return;
        }
        const double abs_det = std::fabs(B.determinant());
        if (abs_det < min_det) {
            return;
        }
        // Direct-cell volume V = 1/|det B|.
        const double direct_vol = 1.0 / abs_det;
        // Protein-ish volume window; blocks noise superlattices and tiny cells.
        if (direct_vol > 2.0e7 || direct_vol < 5.0e3) {
            return;
        }
        const double basis_len = b0.norm() + b1.norm() + b2.norm();
        const Eigen::Matrix3d Binv = lu.inverse();
        const double tol = static_cast<double>(settings.index_tol);
        int explained = 0;
        double index_mag_sum = 0.0;
        int index_n = 0;
        for (const auto& cl : clusters) {
            const Eigen::Vector3d hkl = Binv * cl.center;
            double max_dev = 0.0;
            double mag = 0.0;
            for (int k = 0; k < 3; ++k) {
                const double nearest = std::round(hkl[k]);
                max_dev = std::max(max_dev, std::fabs(hkl[k] - nearest));
                mag += std::fabs(nearest);
            }
            if (max_dev <= tol) {
                explained += cl.count;
                index_mag_sum += mag * static_cast<double>(cl.count);
                index_n += cl.count;
            }
        }
        if (explained < 10 || total_weight <= 0) {
            return;
        }
        const double frac =
            static_cast<double>(explained) / static_cast<double>(total_weight);
        // Need a substantial fraction of cluster weight on-lattice.
        if (frac < 0.35) {
            return;
        }
        const double mean_idx =
            (index_n > 0) ? (index_mag_sum / static_cast<double>(index_n))
                          : 99.0;
        // Rank: (1) explained fraction, (2) shorter reciprocal basis (primitive
        // cell), (3) larger direct volume vs sublattice, (4) smaller indices.
        bool take = false;
        if (frac > best_explained_frac + 0.02) {
            take = true;
        } else if (frac + 0.02 >= best_explained_frac) {
            if (basis_len < best_basis_len * 0.92) {
                take = true;
            } else if (basis_len <= best_basis_len * 1.08) {
                if (direct_vol > best_vol * 1.12) {
                    take = true;
                } else if (
                    direct_vol + 1.0 >= best_vol * 0.88 &&
                    mean_idx < best_mean_idx - 0.05
                ) {
                    take = true;
                }
            }
        }
        if (take) {
            best_explained_frac = frac;
            best_vol = direct_vol;
            best_basis_len = basis_len;
            best_mean_idx = mean_idx;
            best_explained = explained;
            best_om.row(0) = b0;
            best_om.row(1) = b1;
            best_om.row(2) = b2;
        }
    };

    // Prefer short linearly-independent vectors among the dense pool (Kabsch:
    // basis of small magnitude). Try shortest trio first, then all triplets.
    {
        std::vector<std::size_t> by_len = pool_idx;
        std::stable_sort(
            by_len.begin(), by_len.end(),
            [&](std::size_t a, std::size_t b) {
                return clusters[a].center.squaredNorm() <
                       clusters[b].center.squaredNorm();
            }
        );
        std::vector<Eigen::Vector3d> short_basis;
        for (const std::size_t idx : by_len) {
            const Eigen::Vector3d& v = clusters[idx].center;
            if (short_basis.empty()) {
                short_basis.push_back(v);
                continue;
            }
            if (short_basis.size() == 1) {
                if (v.cross(short_basis[0]).squaredNorm() >
                    min_det * min_det) {
                    short_basis.push_back(v);
                }
                continue;
            }
            if (nearly_independent(
                    short_basis[0], short_basis[1], v, min_det
                )) {
                short_basis.push_back(v);
                break;
            }
        }
        if (short_basis.size() == 3) {
            score_basis(short_basis[0], short_basis[1], short_basis[2]);
        }
    }

    for (std::size_t ii = 0; ii < pool_idx.size(); ++ii) {
        for (std::size_t jj = ii + 1; jj < pool_idx.size(); ++jj) {
            for (std::size_t kk = jj + 1; kk < pool_idx.size(); ++kk) {
                const auto& b0 = clusters[pool_idx[ii]].center;
                const auto& b1 = clusters[pool_idx[jj]].center;
                const auto& b2 = clusters[pool_idx[kk]].center;
                if (!nearly_independent(b0, b1, b2, min_det)) {
                    continue;
                }
                score_basis(b0, b1, b2);
            }
        }
    }

    if (best_explained_frac < 0.0 || best_explained < 10) {
        return result;
    }

    // Adjust basis length along each axis by integer factors of cluster
    // vectors. Superlattice (too-large direct cell) ⇒ reciprocal too short
    // ⇒ lengthen; sublattice ⇒ shorten.
    {
        Eigen::Vector3d b[3] = {
            best_om.row(0), best_om.row(1), best_om.row(2)
        };
        Eigen::Matrix3d B0;
        B0.col(0) = b[0];
        B0.col(1) = b[1];
        B0.col(2) = b[2];
        const double vol0 = 1.0 / std::max(1e-18, std::fabs(B0.determinant()));

        auto try_replace = [&](int i, const Eigen::Vector3d& neu) -> bool {
            Eigen::Vector3d trial[3] = {b[0], b[1], b[2]};
            trial[i] = neu;
            if (trial[i].dot(b[i]) < 0.0) {
                trial[i] = -trial[i];
            }
            if (!nearly_independent(trial[0], trial[1], trial[2], min_det)) {
                return false;
            }
            Eigen::Matrix3d Bt;
            Bt.col(0) = trial[0];
            Bt.col(1) = trial[1];
            Bt.col(2) = trial[2];
            const double vol = 1.0 / std::fabs(Bt.determinant());
            if (vol > 2.0e7 || vol < 5.0e3) {
                return false;
            }
            // Do not collapse the cell volume by >2.5× (lengthen overshoot).
            if (vol < vol0 / 2.5) {
                return false;
            }
            Eigen::FullPivLU<Eigen::Matrix3d> lu(Bt);
            if (!lu.isInvertible()) {
                return false;
            }
            const Eigen::Matrix3d Binv = lu.inverse();
            int explained = 0;
            for (const auto& c2 : clusters) {
                const Eigen::Vector3d hkl = Binv * c2.center;
                double max_dev = 0.0;
                for (int k = 0; k < 3; ++k) {
                    max_dev = std::max(
                        max_dev, std::fabs(hkl[k] - std::round(hkl[k]))
                    );
                }
                if (max_dev <= static_cast<double>(settings.index_tol)) {
                    explained += c2.count;
                }
            }
            if (explained + 5 < best_explained) {
                return false;
            }
            b[0] = trial[0];
            b[1] = trial[1];
            b[2] = trial[2];
            best_explained = std::max(best_explained, explained);
            best_om.row(0) = b[0];
            best_om.row(1) = b[1];
            best_om.row(2) = b[2];
            return true;
        };

        bool changed = true;
        int guard = 0;
        while (changed && guard++ < 10) {
            changed = false;
            for (int i = 0; i < 3; ++i) {
                const double bi_len = b[i].norm();
                if (!(bi_len > 1e-12)) {
                    continue;
                }
                const Eigen::Vector3d bi_u = b[i] / bi_len;
                for (const auto& cl : clusters_all) {
                    const double cl_len = cl.center.norm();
                    if (cl_len < len_lo || cl_len > len_hi) {
                        continue;
                    }
                    const Eigen::Vector3d cu = cl.center / cl_len;
                    if (std::fabs(bi_u.dot(cu)) < 0.97) {
                        continue;
                    }
                    // Lengthen: cluster ≈ n * b_i (fix too-large direct edge).
                    {
                        const double ratio = cl_len / bi_len;
                        const int n = static_cast<int>(std::lround(ratio));
                        if (n >= 2 && n <= 8 &&
                            std::fabs(ratio - static_cast<double>(n)) <= 0.12) {
                            if (try_replace(i, cl.center)) {
                                changed = true;
                                break;
                            }
                        }
                    }
                    // Shorten: b_i ≈ n * cluster (fix too-small direct edge).
                    {
                        const double ratio = bi_len / cl_len;
                        const int n = static_cast<int>(std::lround(ratio));
                        if (n >= 2 && n <= 8 &&
                            std::fabs(ratio - static_cast<double>(n)) <= 0.12) {
                            if (try_replace(i, cl.center)) {
                                changed = true;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    result.reciprocal_om = best_om;
    result.cell = cell_from_reciprocal_basis(best_om);
    result.n_explained = best_explained;
    result.ok = result.cell.is_valid();
    return result;
}

}  // namespace mosflm::indexing
