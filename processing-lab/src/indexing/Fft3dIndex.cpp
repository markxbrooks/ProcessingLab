#include "mosflm/indexing/Fft3dIndex.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

namespace mosflm::indexing {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kRad2Deg = 180.0 / kPi;

using Cplx = std::complex<double>;

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

int next_pow2_clamp(int n, int lo, int hi) {
    int p = lo;
    while (p < n && p < hi) {
        p *= 2;
    }
    return std::max(lo, std::min(p, hi));
}

void fft1d(std::vector<Cplx>& a, bool inverse) {
    const int n = static_cast<int>(a.size());
    if (n <= 1) {
        return;
    }
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[static_cast<std::size_t>(i)],
                      a[static_cast<std::size_t>(j)]);
        }
    }
    for (int len = 2; len <= n; len <<= 1) {
        const double ang =
            (inverse ? 2.0 : -2.0) * kPi / static_cast<double>(len);
        const Cplx wlen(std::cos(ang), std::sin(ang));
        for (int i = 0; i < n; i += len) {
            Cplx w(1.0, 0.0);
            for (int j = 0; j < len / 2; ++j) {
                const Cplx u = a[static_cast<std::size_t>(i + j)];
                const Cplx v =
                    a[static_cast<std::size_t>(i + j + len / 2)] * w;
                a[static_cast<std::size_t>(i + j)] = u + v;
                a[static_cast<std::size_t>(i + j + len / 2)] = u - v;
                w *= wlen;
            }
        }
    }
    if (inverse) {
        const double inv_n = 1.0 / static_cast<double>(n);
        for (auto& z : a) {
            z *= inv_n;
        }
    }
}

void fft3d(std::vector<Cplx>& data, int n, bool inverse) {
    std::vector<Cplx> line(static_cast<std::size_t>(n));
    for (int y = 0; y < n; ++y) {
        for (int z = 0; z < n; ++z) {
            for (int x = 0; x < n; ++x) {
                line[static_cast<std::size_t>(x)] =
                    data[static_cast<std::size_t>(x + n * (y + n * z))];
            }
            fft1d(line, inverse);
            for (int x = 0; x < n; ++x) {
                data[static_cast<std::size_t>(x + n * (y + n * z))] =
                    line[static_cast<std::size_t>(x)];
            }
        }
    }
    for (int x = 0; x < n; ++x) {
        for (int z = 0; z < n; ++z) {
            for (int y = 0; y < n; ++y) {
                line[static_cast<std::size_t>(y)] =
                    data[static_cast<std::size_t>(x + n * (y + n * z))];
            }
            fft1d(line, inverse);
            for (int y = 0; y < n; ++y) {
                data[static_cast<std::size_t>(x + n * (y + n * z))] =
                    line[static_cast<std::size_t>(y)];
            }
        }
    }
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < n; ++y) {
            for (int z = 0; z < n; ++z) {
                line[static_cast<std::size_t>(z)] =
                    data[static_cast<std::size_t>(x + n * (y + n * z))];
            }
            fft1d(line, inverse);
            for (int z = 0; z < n; ++z) {
                data[static_cast<std::size_t>(x + n * (y + n * z))] =
                    line[static_cast<std::size_t>(z)];
            }
        }
    }
}

struct Peak {
    Eigen::Vector3d r{Eigen::Vector3d::Zero()};
    double mag{0.0};
};

int score_basis(
    const Eigen::Matrix3d& om,
    const std::vector<Eigen::Vector3d>& rlps,
    double tol
) {
    Eigen::Matrix3d A = om.transpose();
    Eigen::FullPivLU<Eigen::Matrix3d> lu(A);
    if (!lu.isInvertible()) {
        return 0;
    }
    const Eigen::Matrix3d Ainv = lu.inverse();
    int n_ok = 0;
    for (const auto& r : rlps) {
        const Eigen::Vector3d hkl = Ainv * r;
        double max_dev = 0.0;
        for (int j = 0; j < 3; ++j) {
            max_dev =
                std::max(max_dev, std::fabs(hkl[j] - std::round(hkl[j])));
        }
        if (max_dev <= tol) {
            ++n_ok;
        }
    }
    return n_ok;
}

Eigen::Matrix3d reciprocal_from_direct(
    const Eigen::Vector3d& a,
    const Eigen::Vector3d& b,
    const Eigen::Vector3d& c
) {
    const double vol = a.dot(b.cross(c));
    Eigen::Matrix3d om = Eigen::Matrix3d::Zero();
    if (!(std::fabs(vol) > 1e-18)) {
        return om;
    }
    const double inv = 1.0 / vol;
    om.row(0) = b.cross(c) * inv;
    om.row(1) = c.cross(a) * inv;
    om.row(2) = a.cross(b) * inv;
    if (vol < 0.0) {
        om.row(2) = -om.row(2);
    }
    return om;
}

bool nearly_independent(
    const Eigen::Vector3d& a,
    const Eigen::Vector3d& b,
    const Eigen::Vector3d& c,
    double min_abs_det
) {
    return std::fabs(a.dot(b.cross(c))) >= min_abs_det;
}

void sort_edges_len(Eigen::Vector3d& a, Eigen::Vector3d& b, Eigen::Vector3d& c) {
    if (a.norm() > b.norm()) {
        std::swap(a, b);
    }
    if (b.norm() > c.norm()) {
        std::swap(b, c);
    }
    if (a.norm() > b.norm()) {
        std::swap(a, b);
    }
}

bool cells_similar(const UnitCell& u, const UnitCell& v, double rel = 0.05) {
    if (!u.is_valid() || !v.is_valid()) {
        return false;
    }
    double ea[3] = {u.a, u.b, u.c};
    double eb[3] = {v.a, v.b, v.c};
    std::sort(ea, ea + 3);
    std::sort(eb, eb + 3);
    for (int i = 0; i < 3; ++i) {
        const double den = std::max(1e-6, std::fabs(eb[i]));
        if (std::fabs(ea[i] - eb[i]) > rel * den) {
            return false;
        }
    }
    return true;
}

struct Cand {
    Eigen::Matrix3d om{Eigen::Matrix3d::Identity()};
    UnitCell cell{};
    Eigen::Vector3d a{Eigen::Vector3d::Zero()};
    Eigen::Vector3d b{Eigen::Vector3d::Zero()};
    Eigen::Vector3d c{Eigen::Vector3d::Zero()};
    int n_idx{0};
    double vol{0.0};
    double S{0.0};
};

/// DIALS ModelRankFilter::filter_doubled_cell — reject supercells of peers.
std::vector<Cand> filter_doubled_cell(const std::vector<Cand>& solutions) {
    static const int factors[][3] = {
        {2, 1, 1}, {1, 2, 1}, {1, 1, 2}, {2, 2, 1},
        {2, 1, 2}, {1, 2, 2}, {2, 2, 2}
    };
    std::vector<Cand> accepted;
    accepted.reserve(solutions.size());
    for (const auto& s1 : solutions) {
        bool doubled = false;
        for (const auto& f : factors) {
            if (doubled) {
                break;
            }
            UnitCell halved;
            halved.a = s1.cell.a / static_cast<double>(f[0]);
            halved.b = s1.cell.b / static_cast<double>(f[1]);
            halved.c = s1.cell.c / static_cast<double>(f[2]);
            halved.alpha_deg = s1.cell.alpha_deg;
            halved.beta_deg = s1.cell.beta_deg;
            halved.gamma_deg = s1.cell.gamma_deg;
            for (const auto& s2 : solutions) {
                if (&s1 == &s2) {
                    continue;
                }
                if (!cells_similar(halved, s2.cell, 0.05)) {
                    continue;
                }
                if (s1.n_idx < static_cast<int>(1.1 * s2.n_idx + 0.5)) {
                    doubled = true;
                    break;
                }
            }
        }
        if (!doubled) {
            accepted.push_back(s1);
        }
    }
    return accepted;
}

/// DIALS ModelRankFilter, adapted without refinement likelihood:
/// strong n_indexed gate first (avoids tiny junk cells), then doubled-cell,
/// then volume ≤1.25× min among survivors.
std::vector<Cand> model_rank_filter(std::vector<Cand> solutions) {
    if (solutions.empty()) {
        return solutions;
    }
    int max_n = 0;
    for (const auto& s : solutions) {
        max_n = std::max(max_n, s.n_idx);
    }
    // Stricter than DIALS 5% pre-filter: we lack model_likelihood.
    const int pre_cut = std::max(1, static_cast<int>(0.50 * max_n));
    {
        std::vector<Cand> keep;
        for (const auto& s : solutions) {
            if (s.n_idx >= pre_cut) {
                keep.push_back(s);
            }
        }
        solutions = std::move(keep);
    }
    if (solutions.empty()) {
        return solutions;
    }
    max_n = 0;
    for (const auto& s : solutions) {
        max_n = std::max(max_n, s.n_idx);
    }
    const int n_cut = std::max(1, static_cast<int>(0.9 * max_n));
    {
        std::vector<Cand> keep;
        for (const auto& s : solutions) {
            if (s.n_idx >= n_cut) {
                keep.push_back(s);
            }
        }
        solutions = std::move(keep);
    }
    if (solutions.empty()) {
        return solutions;
    }
    solutions = filter_doubled_cell(solutions);
    if (solutions.empty()) {
        return solutions;
    }
    double min_vol = solutions[0].vol;
    for (const auto& s : solutions) {
        min_vol = std::min(min_vol, s.vol);
    }
    {
        std::vector<Cand> keep;
        for (const auto& s : solutions) {
            if (s.vol < 1.25 * min_vol) {
                keep.push_back(s);
            }
        }
        solutions = std::move(keep);
    }
    return solutions;
}

Cand pick_best_cand(const std::vector<Cand>& solutions) {
    auto edge_sum = [](const UnitCell& u) { return u.a + u.b + u.c; };
    Cand best;
    best.n_idx = -1;
    best.vol = 1e300;
    best.S = -1e300;
    for (const auto& s : solutions) {
        if (s.n_idx > best.n_idx) {
            best = s;
            continue;
        }
        if (s.n_idx != best.n_idx) {
            continue;
        }
        // Relative volume tie (FP noise on equivalent cells).
        if (s.vol < best.vol * 0.99) {
            best = s;
            continue;
        }
        if (s.vol <= best.vol * 1.01 &&
            edge_sum(s.cell) < edge_sum(best.cell) - 1e-3) {
            best = s;
        }
    }
    return best;
}

}  // namespace

Fft3dResult find_fft3d_basis(
    const std::vector<float>& x,
    const std::vector<float>& y,
    const std::vector<float>& z,
    double /*wavelength_A*/,
    const Fft3dSettings& settings,
    double dmax_A
) {
    Fft3dResult result;
    const std::size_t n_in = x.size();
    if (n_in < 20 || y.size() != n_in || z.size() != n_in) {
        return result;
    }

    // Real-space search radius (Å) — fixes FFT binning to preserve cell edges.
    double b_max = static_cast<double>(settings.b_max_A);
    if (!(b_max > 10.0)) {
        b_max = (dmax_A > 20.0) ? dmax_A : 200.0;
    }
    const double b_min = std::max(15.0, static_cast<double>(settings.b_min_A));

    // Larger grid when searching long axes (thaumatin c~150 needs ≥128).
    int grid_req = settings.grid_n;
    if (b_max > 100.0 && grid_req < 128) {
        grid_req = 128;
    }
    const int n = next_pow2_clamp(grid_req, 32, 128);
    result.grid_n = n;

    // NEVER enlarge the reciprocal window to fit high-|q| spots — that shrinks
    // max |direct| to ≪ dmax and FFT latches onto ~10 Å junk cells.
    // Δq = 1/(2·b_max); q_half = n·Δq/2 = n/(4·b_max); max |r| = b_max.
    const double dq = 1.0 / (2.0 * b_max);
    const double q_half = 0.5 * static_cast<double>(n) * dq;
    const double inv_voxel = 1.0 / dq;  // = 2·b_max

    const std::size_t n_cap =
        std::min(n_in, std::max<std::size_t>(20, settings.max_spots));
    std::vector<Eigen::Vector3d> rlps_all;
    std::vector<Eigen::Vector3d> rlps_fft;
    rlps_all.reserve(n_cap);
    rlps_fft.reserve(n_cap);
    for (std::size_t i = 0; i < n_cap; ++i) {
        Eigen::Vector3d r(
            static_cast<double>(x[i]),
            static_cast<double>(y[i]),
            static_cast<double>(z[i])
        );
        const double q = r.norm();
        if (q < 1e-8) {
            continue;
        }
        rlps_all.push_back(r);
        // Indexing FFT uses low-|q| RLPs that fit the fixed window.
        if (std::fabs(r.x()) < q_half && std::fabs(r.y()) < q_half &&
            std::fabs(r.z()) < q_half) {
            rlps_fft.push_back(r);
        }
    }
    if (rlps_fft.size() < 20u) {
        return result;
    }
    result.n_spots_used = rlps_fft.size();

    std::vector<Cplx> grid(
        static_cast<std::size_t>(n) * static_cast<std::size_t>(n) *
            static_cast<std::size_t>(n),
        Cplx(0.0, 0.0)
    );
    for (const auto& r : rlps_fft) {
        auto bin = [&](double qi) -> int {
            const double t = (qi + q_half) * inv_voxel;
            int i = static_cast<int>(std::floor(t));
            if (i < 0) {
                i = 0;
            }
            if (i >= n) {
                i = n - 1;
            }
            return i;
        };
        const int ix = bin(r.x());
        const int iy = bin(r.y());
        const int iz = bin(r.z());
        grid[static_cast<std::size_t>(ix + n * (iy + n * iz))] +=
            Cplx(1.0, 0.0);
    }

    fft3d(grid, n, /*inverse=*/false);

    // Real-space: bin k → k·Δr with Δr = 1/(n·Δq) = 2·b_max/n = 1/(2·q_half)·(n/n)
    // wrap(k) / (n·Δq) = wrap(k) · (2·b_max/n) … equivalently wrap/(2·q_half)* (q_half/(n*dq/2))
    // Since q_half = n*dq/2, 2*q_half = n*dq, so wrap/(2*q_half) = wrap/(n*dq) = wrap*Δr. OK.
    auto real_vec = [&](int ix, int iy, int iz) -> Eigen::Vector3d {
        auto wrap = [&](int i) {
            if (i > n / 2) {
                i -= n;
            }
            return i;
        };
        const double sx = static_cast<double>(wrap(ix)) / (2.0 * q_half);
        const double sy = static_cast<double>(wrap(iy)) / (2.0 * q_half);
        const double sz = static_cast<double>(wrap(iz)) / (2.0 * q_half);
        return Eigen::Vector3d(sx, sy, sz);
    };

    std::vector<Peak> peaks;
    peaks.reserve(settings.n_peaks * 4);
    const double origin_cut = b_min * 0.5;
    for (int iz = 0; iz < n; ++iz) {
        for (int iy = 0; iy < n; ++iy) {
            for (int ix = 0; ix < n; ++ix) {
                const double mag = std::abs(
                    grid[static_cast<std::size_t>(ix + n * (iy + n * iz))]
                );
                if (mag < 1e-6) {
                    continue;
                }
                const Eigen::Vector3d rv = real_vec(ix, iy, iz);
                const double len = rv.norm();
                if (len < origin_cut || len > b_max * 1.05) {
                    continue;
                }
                bool is_max = true;
                for (int dz = -1; dz <= 1 && is_max; ++dz) {
                    for (int dy = -1; dy <= 1 && is_max; ++dy) {
                        for (int dx = -1; dx <= 1; ++dx) {
                            if (dx == 0 && dy == 0 && dz == 0) {
                                continue;
                            }
                            const int jx = (ix + dx + n) % n;
                            const int jy = (iy + dy + n) % n;
                            const int jz = (iz + dz + n) % n;
                            const double m2 = std::abs(
                                grid[static_cast<std::size_t>(
                                    jx + n * (jy + n * jz)
                                )]
                            );
                            if (m2 > mag + 1e-12) {
                                is_max = false;
                                break;
                            }
                        }
                    }
                }
                if (!is_max) {
                    continue;
                }
                Peak p;
                p.r = fold_hemisphere(rv);
                p.mag = mag;
                bool dup = false;
                for (const auto& q : peaks) {
                    if ((p.r - q.r).norm() < b_min * 0.2) {
                        dup = true;
                        break;
                    }
                }
                if (!dup) {
                    peaks.push_back(p);
                }
            }
        }
    }

    std::sort(peaks.begin(), peaks.end(), [](const Peak& a, const Peak& b) {
        return a.mag > b.mag;
    });
    if (peaks.size() > settings.n_peaks) {
        peaks.resize(settings.n_peaks);
    }
    result.n_peaks = peaks.size();
    if (peaks.size() < 3u) {
        return result;
    }
    if (std::getenv("FFT3D_DEBUG_PEAKS") != nullptr) {
        std::cout << "  FFT3D peaks (" << peaks.size() << "):";
        for (std::size_t i = 0; i < peaks.size() && i < 20; ++i) {
            std::cout << " " << std::fixed << peaks[i].r.norm()
                      << "(" << peaks[i].mag << ")";
        }
        std::cout << '\n';
    }

    // Inject (u±v)/2 for near-equal angled pairs — recovers primitive edges
    // when FFT prefers face diagonals (~58√2 → ~58 for thaumatin).
    auto add_peak = [&](Eigen::Vector3d r, double mag) {
        r = fold_hemisphere(std::move(r));
        const double len = r.norm();
        if (len < b_min || len > b_max * 1.05) {
            return;
        }
        for (const auto& q : peaks) {
            if ((r - q.r).norm() < b_min * 0.15) {
                return;
            }
        }
        Peak p;
        p.r = r;
        p.mag = mag;
        peaks.push_back(p);
    };
    const std::size_t n_fft_peaks = peaks.size();
    for (std::size_t i = 0; i < n_fft_peaks; ++i) {
        for (std::size_t j = i + 1; j < n_fft_peaks; ++j) {
            const Eigen::Vector3d& u = peaks[i].r;
            const Eigen::Vector3d& v = peaks[j].r;
            const double lu = u.norm();
            const double lv = v.norm();
            if (lu < 1e-6 || lv < 1e-6) {
                continue;
            }
            const double ratio = (lu > lv) ? (lu / lv) : (lv / lu);
            if (ratio > 1.15) {
                continue;
            }
            const double cang =
                std::max(-1.0, std::min(1.0, u.dot(v) / (lu * lv)));
            const double ang = std::acos(cang) * kRad2Deg;
            if (ang < 60.0 || ang > 120.0) {
                continue;
            }
            const double mag = 0.5 * (peaks[i].mag + peaks[j].mag);
            add_peak(0.5 * (u + v), mag);
            add_peak(0.5 * (u - v), mag);
        }
    }

    // Mix strongest FFT peaks with shortest (primitive) candidates.
    std::vector<Peak> ranked = peaks;
    std::vector<Peak> by_len = peaks;
    std::sort(by_len.begin(), by_len.end(), [](const Peak& a, const Peak& b) {
        return a.r.norm() < b.r.norm();
    });
    {
        std::vector<Peak> mixed;
        mixed.reserve(40);
        for (std::size_t i = 0; i < ranked.size() && mixed.size() < 20; ++i) {
            mixed.push_back(ranked[i]);
        }
        for (std::size_t i = 0; i < by_len.size() && mixed.size() < 36; ++i) {
            bool seen = false;
            for (const auto& m : mixed) {
                if ((m.r - by_len[i].r).norm() < 1.0) {
                    seen = true;
                    break;
                }
            }
            if (!seen) {
                mixed.push_back(by_len[i]);
            }
        }
        ranked = std::move(mixed);
    }

    // Slightly tighter than CONTROL default — supercells inflate loose fits.
    const double tol = std::min(0.22, static_cast<double>(settings.index_tol));
    const double min_vol = std::max(5000.0, 0.02 * b_max * b_max * b_max);
    const double V0 = b_max * b_max * b_max;
    constexpr double kVolLambda = 100.0;
    auto vol_score = [&](int n_idx, double vol) -> double {
        const double vv = std::max(vol, 1.0);
        return static_cast<double>(n_idx) - kVolLambda * std::log(vv / V0);
    };

    std::vector<Cand> solutions;
    solutions.reserve(48);

    auto try_add = [&](Eigen::Vector3d a, Eigen::Vector3d b,
                       Eigen::Vector3d c) {
        sort_edges_len(a, b, c);
        if (a.norm() < b_min) {
            return;
        }
        const double triple = std::fabs(a.dot(b.cross(c)));
        if (triple < min_vol) {
            return;
        }
        // Reject near-coplanar bases (flat parallelepiped, tiny triple vs edges).
        const double edge_prod = a.norm() * b.norm() * c.norm();
        if (edge_prod > 1e-6 && triple < 0.20 * edge_prod) {
            return;
        }
        const Eigen::Matrix3d om = reciprocal_from_direct(a, b, c);
        if (om.squaredNorm() < 1e-18) {
            return;
        }
        const int sc = score_basis(om, rlps_fft, tol);
        // Soft floor: synthetic lattices clear easily; real spots are noisier.
        const int min_sc = std::max(
            10, static_cast<int>(0.12 * static_cast<double>(rlps_fft.size()))
        );
        if (sc < min_sc) {
            return;
        }
        Cand cand;
        cand.om = om;
        cand.cell = cell_from_reciprocal_basis(om);
        if (!cand.cell.is_valid()) {
            return;
        }
        if (cand.cell.a < b_min * 0.8 || cand.cell.b < b_min * 0.8 ||
            cand.cell.c < b_min * 0.8) {
            return;
        }
        // Rank by crystallographic cell volume (edges×angles), not a possibly
        // flat triple product that under-reports supercell volume.
        const double vol = std::fabs(cand.cell.volume());
        if (vol < min_vol) {
            return;
        }
        cand.a = a;
        cand.b = b;
        cand.c = c;
        cand.n_idx = sc;
        cand.vol = vol;
        cand.S = vol_score(sc, vol);
        for (const auto& s : solutions) {
            if (cells_similar(s.cell, cand.cell, 0.03) &&
                cand.n_idx <= s.n_idx && cand.vol >= 0.98 * s.vol) {
                return;
            }
        }
        solutions.push_back(std::move(cand));
    };

    // DIALS: three shortest non-collinear peaks with appreciable weight.
    // Restrict to the stronger half of the peak list, then sort by length.
    // Independence floor is softer than min_vol — try_add still enforces volume.
    const double min_cross2 = b_min * b_min * 0.05;
    const double min_det = std::max(500.0, 0.01 * min_vol);
    {
        std::vector<Peak> strong = ranked;
        if (strong.size() > 8) {
            const std::size_t keep =
                std::max<std::size_t>(8, strong.size() / 2);
            strong.resize(keep);
        }
        std::sort(
            strong.begin(), strong.end(),
            [](const Peak& x, const Peak& y) {
                return x.r.norm() < y.r.norm();
            }
        );
        std::vector<Eigen::Vector3d> short_basis;
        for (const auto& p : strong) {
            const Eigen::Vector3d& v = p.r;
            if (v.norm() < b_min) {
                continue;
            }
            if (short_basis.empty()) {
                short_basis.push_back(v);
                continue;
            }
            if (short_basis.size() == 1) {
                if (v.cross(short_basis[0]).squaredNorm() > min_cross2) {
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
        if (std::getenv("FFT3D_DEBUG_PEAKS") != nullptr) {
            std::cout << "  FFT3D short_basis:";
            for (const auto& v : short_basis) {
                std::cout << " " << v.norm();
            }
            std::cout << '\n';
        }
        if (short_basis.size() == 3) {
            try_add(short_basis[0], short_basis[1], short_basis[2]);
        }
    }

    // Shortest-first triplets (primitive edges before diagonal supercells),
    // then a few strong-peak combos for diversity.
    auto enum_triplets = [&](const std::vector<Peak>& pool_peaks,
                             std::size_t cap) {
        const std::size_t n_p = std::min(pool_peaks.size(), cap);
        for (std::size_t i = 0; i < n_p; ++i) {
            for (std::size_t j = i + 1; j < n_p; ++j) {
                for (std::size_t k = j + 1; k < n_p; ++k) {
                    try_add(
                        pool_peaks[i].r, pool_peaks[j].r, pool_peaks[k].r
                    );
                }
            }
        }
    };
    enum_triplets(by_len, 24);
    enum_triplets(ranked, 16);

    if (solutions.empty()) {
        return result;
    }

    // Prefer high n_idx, then smaller volume, then shorter edge sum.
    std::sort(solutions.begin(), solutions.end(), [](const Cand& x,
                                                     const Cand& y) {
        if (x.n_idx != y.n_idx) {
            return x.n_idx > y.n_idx;
        }
        if (x.vol < y.vol * 0.99) {
            return true;
        }
        if (y.vol < x.vol * 0.99) {
            return false;
        }
        return (x.cell.a + x.cell.b + x.cell.c) <
               (y.cell.a + y.cell.b + y.cell.c);
    });
    if (solutions.size() > 40) {
        solutions.resize(40);
    }

    if (std::getenv("FFT3D_DEBUG_PEAKS") != nullptr) {
        std::cout << "  FFT3D solutions (" << solutions.size() << "):";
        for (std::size_t i = 0; i < solutions.size() && i < 12; ++i) {
            const auto& s = solutions[i];
            std::cout << " [" << s.cell.a << " " << s.cell.b << " " << s.cell.c
                      << " n=" << s.n_idx << " V=" << s.vol << "]";
        }
        std::cout << '\n';
    }

    solutions = model_rank_filter(std::move(solutions));
    if (solutions.empty()) {
        return result;
    }

    if (std::getenv("FFT3D_DEBUG_PEAKS") != nullptr) {
        std::cout << "  FFT3D filtered (" << solutions.size() << "):";
        for (const auto& s : solutions) {
            std::cout << " [" << s.cell.a << " " << s.cell.b << " " << s.cell.c
                      << " n=" << s.n_idx << " V=" << s.vol << "]";
        }
        std::cout << '\n';
    }

    Cand best = pick_best_cand(solutions);
    int best_nidx = best.n_idx;
    double best_vol = best.vol;
    double best_S = best.S;
    Eigen::Matrix3d best_om = best.om;
    UnitCell best_cell = best.cell;
    Eigen::Vector3d a = best.a;
    Eigen::Vector3d b = best.b;
    Eigen::Vector3d c = best.c;

    // Reduce supercells: basal √2, (vi±vj)/2, axis-halving, and short
    // lattice vectors i·a+j·b+k·c (|i|,|j|,|k|≤2) that still index well.
    {
        const double vol_before = best_vol;
        const UnitCell cell_before = best_cell;

        auto edge_sum = [](const UnitCell& u) {
            return u.a + u.b + u.c;
        };
        auto accept_reduced = [&](int sc2, double vol2,
                                  const UnitCell& cell2) -> bool {
            if (sc2 + 1 < static_cast<int>(0.9 * best_nidx)) {
                return false;
            }
            if (sc2 < best_nidx - 2) {
                return false;
            }
            // Smaller volume (supercell → primitive).
            if (vol2 < 0.92 * best_vol) {
                return true;
            }
            // Unimodular reorientation: same volume, shorter edges (Niggli-ish).
            if (std::fabs(vol2 - best_vol) <= 0.08 * best_vol &&
                edge_sum(cell2) < 0.95 * edge_sum(best_cell)) {
                return true;
            }
            return false;
        };

        auto consider = [&](Eigen::Vector3d t0, Eigen::Vector3d t1,
                            Eigen::Vector3d t2) -> bool {
            sort_edges_len(t0, t1, t2);
            if (t0.norm() < b_min || t1.norm() < b_min ||
                t2.norm() < b_min) {
                return false;
            }
            if (!nearly_independent(t0, t1, t2, min_vol)) {
                return false;
            }
            const double triple2 = std::fabs(t0.dot(t1.cross(t2)));
            const double edge_prod2 = t0.norm() * t1.norm() * t2.norm();
            if (triple2 < min_vol ||
                (edge_prod2 > 1e-6 && triple2 < 0.20 * edge_prod2)) {
                return false;
            }
            const Eigen::Matrix3d om2 =
                reciprocal_from_direct(t0, t1, t2);
            const UnitCell cell2 = cell_from_reciprocal_basis(om2);
            const double vol2 = std::fabs(cell2.volume());
            if (vol2 < min_vol) {
                return false;
            }
            const int sc2 = score_basis(om2, rlps_fft, tol);
            if (!accept_reduced(sc2, vol2, cell2)) {
                return false;
            }
            a = t0;
            b = t1;
            c = t2;
            best_om = om2;
            best_cell = cell2;
            best_nidx = sc2;
            best_vol = vol2;
            best_S = vol_score(sc2, vol2);
            return true;
        };

        for (int pass = 0; pass < 16; ++pass) {
            Eigen::Vector3d vecs[3] = {a, b, c};
            bool changed = false;

            // Explicit basal √2: {u,v} → {(u+v)/2,(u−v)/2}.
            for (int i = 0; i < 3; ++i) {
                for (int j = i + 1; j < 3; ++j) {
                    const double li = vecs[i].norm();
                    const double lj = vecs[j].norm();
                    if (li < 1e-6 || lj < 1e-6) {
                        continue;
                    }
                    const double ratio =
                        (li > lj) ? (li / lj) : (lj / li);
                    if (ratio > 1.15) {
                        continue;
                    }
                    const double cang = std::max(
                        -1.0,
                        std::min(1.0, vecs[i].dot(vecs[j]) / (li * lj))
                    );
                    const double ang = std::acos(cang) * kRad2Deg;
                    if (ang < 60.0 || ang > 120.0) {
                        continue;
                    }
                    const int k = 3 - i - j;
                    Eigen::Vector3d trial[3];
                    trial[i] = 0.5 * (vecs[i] + vecs[j]);
                    trial[j] = 0.5 * (vecs[i] - vecs[j]);
                    trial[k] = vecs[k];
                    if (consider(trial[0], trial[1], trial[2])) {
                        changed = true;
                        vecs[0] = a;
                        vecs[1] = b;
                        vecs[2] = c;
                    }
                }
            }

            // DIALS-style axis halving (doubled / tripled cells).
            for (int m0 = 1; m0 <= 2; ++m0) {
                for (int m1 = 1; m1 <= 2; ++m1) {
                    for (int m2 = 1; m2 <= 2; ++m2) {
                        if (m0 * m1 * m2 == 1) {
                            continue;
                        }
                        Eigen::Vector3d trial[3] = {
                            vecs[0] / static_cast<double>(m0),
                            vecs[1] / static_cast<double>(m1),
                            vecs[2] / static_cast<double>(m2)
                        };
                        if (consider(trial[0], trial[1], trial[2])) {
                            changed = true;
                            vecs[0] = a;
                            vecs[1] = b;
                            vecs[2] = c;
                        }
                    }
                }
            }

            // Replace one vector by a shorter lattice combo i·a+j·b+k·c.
            for (int i = -2; i <= 2; ++i) {
                for (int j = -2; j <= 2; ++j) {
                    for (int k = -2; k <= 2; ++k) {
                        if (i == 0 && j == 0 && k == 0) {
                            continue;
                        }
                        // Skip pure ± basis (no change) and gcd>1 scales
                        // handled by axis-halving.
                        const int abs_sum =
                            std::abs(i) + std::abs(j) + std::abs(k);
                        if (abs_sum == 1) {
                            continue;
                        }
                        const Eigen::Vector3d t =
                            static_cast<double>(i) * vecs[0] +
                            static_cast<double>(j) * vecs[1] +
                            static_cast<double>(k) * vecs[2];
                        const double tlen = t.norm();
                        if (tlen < b_min || tlen > b_max * 1.05) {
                            continue;
                        }
                        for (int axis = 0; axis < 3; ++axis) {
                            if (!(tlen + 0.5 < vecs[axis].norm())) {
                                continue;
                            }
                            Eigen::Vector3d trial[3] = {
                                vecs[0], vecs[1], vecs[2]
                            };
                            trial[axis] = t;
                            if (consider(trial[0], trial[1], trial[2])) {
                                changed = true;
                                vecs[0] = a;
                                vecs[1] = b;
                                vecs[2] = c;
                            }
                        }
                    }
                }
            }

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    if (i == j) {
                        continue;
                    }
                    for (double s : {1.0, -1.0}) {
                        for (double scale : {1.0, 0.5}) {
                            const Eigen::Vector3d t =
                                scale * (vecs[i] + s * vecs[j]);
                            const double tlen = t.norm();
                            if (tlen < b_min || tlen > b_max * 1.05) {
                                continue;
                            }
                            if (scale >= 1.0 &&
                                !(tlen + 0.5 < vecs[i].norm())) {
                                continue;
                            }
                            Eigen::Vector3d trial[3] = {
                                vecs[0], vecs[1], vecs[2]
                            };
                            trial[i] = t;
                            if (consider(trial[0], trial[1], trial[2])) {
                                changed = true;
                                vecs[0] = a;
                                vecs[1] = b;
                                vecs[2] = c;
                            }
                        }
                    }
                }
            }
            if (!changed) {
                break;
            }
        }
        if (best_vol < 0.90 * vol_before) {
            std::cout << "  FFT3D: reduced " << std::fixed
                      << cell_before.a << " " << cell_before.b << " "
                      << cell_before.c << " → " << best_cell.a << " "
                      << best_cell.b << " " << best_cell.c << '\n';
        }
    }

    // Report indexed count on the full capped spot list.
    result.reciprocal_om = best_om;
    result.cell = best_cell;
    result.n_indexed = score_basis(best_om, rlps_all, tol);
    result.ok = true;
    return result;
}

}  // namespace mosflm::indexing
