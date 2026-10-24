#include "mosflm/diffraction/Strategy.hpp"

#include "mosflm/geometry/ReciprocalLattice.hpp"

#include <gemmi/symmetry.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>

namespace mosflm::diffraction {
namespace {

constexpr double kRad = 3.14159265358979323846 / 180.0;

std::int64_t hkl_key(int h, int k, int l) {
    return (static_cast<std::int64_t>(h + 1024) << 22) | (static_cast<std::int64_t>(k + 1024) << 11) |
           static_cast<std::int64_t>(l + 1024);
}

/// A unique reflection and the φ bins it is recorded in (each mate).
struct Unique {
    double d_A{0.0};
    bool centric{false};
    std::vector<int> plus;   ///< Bins of I+ (and centric) observations.
    std::vector<int> minus;  ///< Bins of I- observations.
    std::vector<int> any;    ///< Both, sorted.
};

using Intervals = std::vector<std::pair<int, int>>;

/// Window starts ``s`` in ``[0, max_start]`` whose window ``[s, s + width)``
/// holds one of ``bins`` (sorted), as merged intervals.
void starts_covering(const std::vector<int>& bins, int width, int max_start, Intervals& out) {
    out.clear();
    for (const int b : bins) {
        const int lo = std::max(0, b - width + 1);
        const int hi = std::min(b, max_start);
        if (lo > hi) {
            continue;
        }
        if (!out.empty() && lo <= out.back().second + 1) {
            out.back().second = std::max(out.back().second, hi);
        } else {
            out.emplace_back(lo, hi);
        }
    }
}

void intersect(const Intervals& a, const Intervals& b, Intervals& out) {
    out.clear();
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < a.size() && j < b.size()) {
        const int lo = std::max(a[i].first, b[j].first);
        const int hi = std::min(a[i].second, b[j].second);
        if (lo <= hi) {
            out.emplace_back(lo, hi);
        }
        (a[i].second < b[j].second ? i : j)++;
    }
}

}  // namespace

StrategyResult compute_strategy(const PredictionRequest& geometry, const StrategySettings& settings) {
    const auto* sg = geometry.crystal.spacegroup;
    if (sg == nullptr) {
        throw std::invalid_argument("STRATEGY needs a space group");
    }
    const auto& ps = geometry.settings;
    const auto& det = geometry.detector;
    const double lambda = det.wavelength_A;
    if (!(lambda > 0.0) || !(ps.dmin_A > 0.0)) {
        throw std::invalid_argument("STRATEGY needs the wavelength and a resolution limit");
    }
    const double step = std::max(settings.step_deg, 0.1);
    const UnitCell cell{geometry.crystal.cell.a,     geometry.crystal.cell.b,    geometry.crystal.cell.c,
                        geometry.crystal.cell.alpha, geometry.crystal.cell.beta, geometry.crystal.cell.gamma};
    const ReciprocalLattice recip = ReciprocalLattice::from_direct(cell);
    const gemmi::ReciprocalAsu asu(sg);
    const gemmi::GroupOps gops = sg->operations();

    StrategyResult res;
    res.dmin_A = ps.dmin_A;
    res.dmax_A = ps.dmax_A;
    res.osc_deg = settings.osc_deg;

    // Every unique reflection the resolution range allows.
    const int hmax = static_cast<int>(std::ceil(cell.a / ps.dmin_A)) + 1;
    const int kmax = static_cast<int>(std::ceil(cell.b / ps.dmin_A)) + 1;
    const int lmax = static_cast<int>(std::ceil(cell.c / ps.dmin_A)) + 1;
    std::vector<Unique> uniques;
    std::unordered_map<std::int64_t, int> index;
    for (int h = -hmax; h <= hmax; ++h) {
        for (int k = -kmax; k <= kmax; ++k) {
            for (int l = -lmax; l <= lmax; ++l) {
                const gemmi::Op::Miller m{h, k, l};
                if ((h == 0 && k == 0 && l == 0) || !asu.is_in(m) || gops.is_systematically_absent(m)) {
                    continue;
                }
                const double d = recip.d_spacing(h, k, l);
                if (d < ps.dmin_A || d > ps.dmax_A) {
                    continue;
                }
                index.emplace(hkl_key(h, k, l), static_cast<int>(uniques.size()));
                Unique u;
                u.d_A = d;
                u.centric = gops.is_reflection_centric(m);
                uniques.push_back(std::move(u));
            }
        }
    }
    res.unique_possible = static_cast<int>(uniques.size());
    res.acentric_possible = static_cast<int>(
        std::count_if(uniques.begin(), uniques.end(), [](const Unique& u) { return !u.centric; }));
    if (uniques.empty()) {
        return res;
    }
    res.dmax_A = std::max_element(uniques.begin(), uniques.end(), [](const Unique& x, const Unique& y) {
                     return x.d_A < y.d_A;
                 })->d_A;

    // φ domain: bins of ``step`` from ``origin``. A full turn wraps (bins
    // repeated one turn on); a collected sweep doesn't.
    const bool circular = !settings.sweep_deg.has_value();
    const double origin = circular ? 0.0 : settings.sweep_deg->first;
    const double span = circular ? 360.0 : std::min(settings.sweep_deg->second - origin, 360.0);
    const int n_bins = std::max(1, static_cast<int>(std::lround(span / step)));
    const auto bin_of = [&](double phi) -> int {
        double rel = std::fmod(phi - origin, 360.0);
        if (rel < 0.0) {
            rel += 360.0;
        }
        const int b = static_cast<int>(std::floor(rel / step + 1e-9));
        return b < n_bins ? b : -1;
    };

    // Follow every reflection over a full turn (both Ewald crossings).
    const Eigen::Matrix3d& ub = geometry.orientation.ub();
    const Eigen::Vector3d s0(0.0, 0.0, 1.0 / lambda);
    const double margin = ps.apply_reek_detector_margin ? ps.detector_margin_px : 0.0;
    for (int h = -hmax; h <= hmax; ++h) {
        for (int k = -kmax; k <= kmax; ++k) {
            for (int l = -lmax; l <= lmax; ++l) {
                if (h == 0 && k == 0 && l == 0) {
                    continue;
                }
                const gemmi::Op::Miller m{h, k, l};
                const auto [am, isym] = asu.to_asu(m, gops);
                const auto it = index.find(hkl_key(am[0], am[1], am[2]));
                if (it == index.end()) {
                    continue;  // outside the resolution range, or absent
                }
                const Eigen::Vector3d q0 = ub * Eigen::Vector3d(h, k, l);
                const double qn = q0.norm();
                // q_y sinφ + q_z cosφ = −(λ/2)|q|²  (rotation about +X).
                const double amp = std::hypot(q0.y(), q0.z());
                const double rhs = -0.5 * lambda * qn * qn;
                if (amp < 1e-15 || std::fabs(rhs) > amp) {
                    continue;  // never crosses the sphere: the blind region
                }
                const double alpha = std::atan2(q0.y(), q0.z());
                const double delta = std::acos(std::clamp(rhs / amp, -1.0, 1.0));
                for (const double sol : {alpha + delta, alpha - delta}) {
                    const double phi = sol / kRad;
                    const double c = std::cos(sol);
                    const double s = std::sin(sol);
                    const Eigen::Vector3d q(q0.x(), c * q0.y() - s * q0.z(), s * q0.y() + c * q0.z());
                    const Eigen::Vector3d s1 = s0 + q;
                    const Eigen::Vector3d e1 = s1.cross(s0);
                    if (e1.norm() < 1e-12 || std::fabs(e1.x()) / e1.norm() < std::max(ps.min_zeta, 1e-6)) {
                        continue;  // too close to the spindle to measure
                    }
                    double x = 0.0;
                    double y = 0.0;
                    if (!det.lab_to_detector_mm(s1, x, y) ||
                        !spot_on_panel(det, x, y, margin, ps.use_panel_beam, ps.panel_beam_x_mm,
                                       ps.panel_beam_y_mm)) {
                        continue;
                    }
                    const int b = bin_of(phi);
                    if (b < 0) {
                        continue;
                    }
                    Unique& u = uniques[static_cast<std::size_t>(it->second)];
                    // ISYM odd: I+ (or centric), even: I-.
                    (u.centric || isym % 2 == 1 ? u.plus : u.minus).push_back(b);
                }
            }
        }
    }
    for (auto& u : uniques) {
        if (circular) {
            // A window may wrap past 360°: repeat each bin one turn on.
            for (auto* v : {&u.plus, &u.minus}) {
                const std::size_t n = v->size();
                for (std::size_t i = 0; i < n; ++i) {
                    v->push_back((*v)[i] + n_bins);
                }
            }
        }
        std::sort(u.plus.begin(), u.plus.end());
        std::sort(u.minus.begin(), u.minus.end());
        u.any = u.plus;
        u.any.insert(u.any.end(), u.minus.begin(), u.minus.end());
        std::sort(u.any.begin(), u.any.end());
    }

    const double n_unique = static_cast<double>(res.unique_possible);
    const double n_acentric = std::max(1.0, static_cast<double>(res.acentric_possible));
    const auto segment = [&](int start, int width) {
        StrategySegment seg;
        seg.start_deg = origin + start * step;
        seg.rotation_deg = width * step;
        seg.end_deg = seg.start_deg + seg.rotation_deg;
        if (circular) {
            seg.start_deg = std::fmod(seg.start_deg, 360.0);
            seg.end_deg = seg.start_deg + seg.rotation_deg;
        }
        seg.images = static_cast<int>(std::ceil(seg.rotation_deg / std::max(settings.osc_deg, 1e-6) - 1e-6));
        const auto in = [&](const std::vector<int>& bins) {
            return static_cast<int>(std::count_if(bins.begin(), bins.end(),
                                                  [&](int b) { return b >= start && b < start + width; }));
        };
        // Shells of equal reciprocal volume, from the lowest resolution recorded.
        const int ns = std::max(1, settings.n_shells);
        double dlow = 0.0;
        for (const auto& u : uniques) {
            dlow = std::max(dlow, u.d_A);
        }
        const double s_lo = 1.0 / (dlow * dlow * dlow);
        const double s_hi = 1.0 / (ps.dmin_A * ps.dmin_A * ps.dmin_A);
        seg.shells.resize(static_cast<std::size_t>(ns));
        for (int i = 0; i < ns; ++i) {
            auto& sh = seg.shells[static_cast<std::size_t>(i)];
            sh.d_low_A = std::cbrt(1.0 / (s_lo + (s_hi - s_lo) * i / ns));
            sh.d_high_A = std::cbrt(1.0 / (s_lo + (s_hi - s_lo) * (i + 1) / ns));
        }
        std::vector<int> obs(static_cast<std::size_t>(ns), 0);
        std::vector<int> acentric(static_cast<std::size_t>(ns), 0);
        std::vector<int> pairs(static_cast<std::size_t>(ns), 0);
        int total_observed = 0;
        int total_obs = 0;
        int total_pairs = 0;
        for (const auto& u : uniques) {
            const double s3 = 1.0 / (u.d_A * u.d_A * u.d_A);
            const int shell = std::clamp(static_cast<int>((s3 - s_lo) / (s_hi - s_lo) * ns), 0, ns - 1);
            auto& sh = seg.shells[static_cast<std::size_t>(shell)];
            ++sh.possible;
            const int np = in(u.plus);
            const int nm = in(u.minus);
            if (np + nm > 0) {
                ++sh.observed;
                ++total_observed;
                obs[static_cast<std::size_t>(shell)] += np + nm;
                total_obs += np + nm;
            }
            if (!u.centric) {
                ++acentric[static_cast<std::size_t>(shell)];
                if (np > 0 && nm > 0) {
                    ++pairs[static_cast<std::size_t>(shell)];
                    ++total_pairs;
                }
            }
        }
        for (int i = 0; i < ns; ++i) {
            auto& sh = seg.shells[static_cast<std::size_t>(i)];
            const auto k = static_cast<std::size_t>(i);
            sh.completeness = sh.possible > 0 ? static_cast<double>(sh.observed) / sh.possible : 0.0;
            sh.multiplicity = sh.observed > 0 ? static_cast<double>(obs[k]) / sh.observed : 0.0;
            sh.anomalous_completeness = acentric[k] > 0 ? static_cast<double>(pairs[k]) / acentric[k] : 0.0;
        }
        seg.completeness = total_observed / n_unique;
        seg.anomalous_completeness = total_pairs / n_acentric;
        seg.multiplicity = total_observed > 0 ? static_cast<double>(total_obs) / total_observed : 0.0;
        return seg;
    };

    // For each rotation (in bins), the window start recording the most.
    const int max_width = n_bins;
    const int wanted = settings.rotation_deg > 0.0
                           ? std::clamp(static_cast<int>(std::lround(settings.rotation_deg / step)), 1, max_width)
                           : 0;
    std::vector<int> diff_any;
    std::vector<int> diff_anom;
    Intervals a;
    Intervals p;
    Intervals mi;
    Intervals both;
    for (int width = 1; width <= max_width; ++width) {
        const int max_start = circular ? n_bins - 1 : n_bins - width;
        diff_any.assign(static_cast<std::size_t>(max_start + 2), 0);
        diff_anom.assign(static_cast<std::size_t>(max_start + 2), 0);
        for (const auto& u : uniques) {
            starts_covering(u.any, width, max_start, a);
            for (const auto& [lo, hi] : a) {
                ++diff_any[static_cast<std::size_t>(lo)];
                --diff_any[static_cast<std::size_t>(hi + 1)];
            }
            if (!u.centric && !u.plus.empty() && !u.minus.empty()) {
                starts_covering(u.plus, width, max_start, p);
                starts_covering(u.minus, width, max_start, mi);
                intersect(p, mi, both);
                for (const auto& [lo, hi] : both) {
                    ++diff_anom[static_cast<std::size_t>(lo)];
                    --diff_anom[static_cast<std::size_t>(hi + 1)];
                }
            }
        }
        // The best start for the quantity the target is about; ties: earliest.
        int best = 0;
        int best_any = -1;
        int best_anom = -1;
        int run_any = 0;
        int run_anom = 0;
        for (int s = 0; s <= max_start; ++s) {
            run_any += diff_any[static_cast<std::size_t>(s)];
            run_anom += diff_anom[static_cast<std::size_t>(s)];
            const bool better = settings.anomalous
                                    ? std::tie(run_anom, run_any) > std::tie(best_anom, best_any)
                                    : std::tie(run_any, run_anom) > std::tie(best_any, best_anom);
            if (better) {
                best = s;
                best_any = run_any;
                best_anom = run_anom;
            }
        }
        StrategyPoint pt;
        pt.rotation_deg = width * step;
        pt.start_deg = std::fmod(origin + best * step, 360.0);
        pt.completeness = best_any / n_unique;
        pt.anomalous_completeness = best_anom / n_acentric;
        res.curve.push_back(pt);
        const double reached = settings.anomalous ? pt.anomalous_completeness : pt.completeness;
        if (!res.target && reached >= settings.target_completeness - 1e-12) {
            res.target = segment(best, width);
        }
        if (width == wanted) {
            res.rotation = segment(best, width);
        }
        // Nothing more to gain once every reflection that can be recorded is.
        if (res.target && (wanted == 0 || width >= wanted) && pt.completeness >= 1.0 - 1e-12 &&
            (!settings.anomalous || pt.anomalous_completeness >= 1.0 - 1e-12)) {
            break;
        }
    }
    if (!res.curve.empty()) {
        res.max_completeness = res.curve.back().completeness;
        res.max_anomalous_completeness = res.curve.back().anomalous_completeness;
        for (const auto& pt : res.curve) {
            res.max_completeness = std::max(res.max_completeness, pt.completeness);
            res.max_anomalous_completeness = std::max(res.max_anomalous_completeness, pt.anomalous_completeness);
        }
    }
    return res;
}

}  // namespace mosflm::diffraction
