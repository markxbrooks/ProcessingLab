#include "mosflm/indexing/AutoIndexer.hpp"

#include "mosflm/indexing/BravaisSelect.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

extern "C" {
#include "util.h"
#include "index_def.h"
#include "cryst_matrix.h"
}

namespace mosflm::indexing {
namespace {

extern "C" void index_fortran_(
    int* numvect,
    int* numvectref,
    int* maxcyc,
    int* inicyc,
    float* convrg,
    int* numormat,
    float* stepscal,
    int* numrlp,
    int* numsptdps,
    int* npkfftmin,
    int* nmultfft,
    int* idumpth,
    float* theta1d,
    float* theta2d,
    float* sinthmin,
    float* volmin,
    float* hkldev,
    float* hklmin,
    int dumpfft[],
    int* dumpvec,
    int* dumpref,
    int* dumpgen,
    int* dumpch,
    int* number_xyz,
    int* nlattice,
    float* x,
    float* y,
    float* z,
    int* usespot,
    float* d_max,
    float* d_min,
    int* index_mode,
    float* a,
    float* wave,
    int* ierrflg,
    float* tarcell,
    int* centering
);

UnitCell cell_from_dps(const struct cell_struct& c) {
    UnitCell out;
    out.a = c.a;
    out.b = c.b;
    out.c = c.c;
    constexpr double kRad = 180.0 / 3.14159265358979323846;
    out.alpha_deg = c.alpha * kRad;
    out.beta_deg = c.beta * kRad;
    out.gamma_deg = c.gamma * kRad;
    return out;
}

Eigen::Matrix3d amat_from_dps_om(const struct matrix_3x3& om, double wave) {
    // ind_fortran.c: MOSFLM A = λ · permutation of DPS OM rows.
    Eigen::Matrix3d A;
    A(0, 0) = wave * om.r[0].z;
    A(0, 1) = wave * om.r[0].x;
    A(0, 2) = wave * om.r[0].y;
    A(1, 0) = wave * om.r[1].z;
    A(1, 1) = wave * om.r[1].x;
    A(1, 2) = wave * om.r[1].y;
    A(2, 0) = wave * om.r[2].z;
    A(2, 1) = wave * om.r[2].x;
    A(2, 2) = wave * om.r[2].y;
    return A;
}

struct matrix_3x3 eigen_to_dps(const Eigen::Matrix3d& M) {
    struct matrix_3x3 om{};
    for (int i = 0; i < 3; ++i) {
        om.r[i].x = static_cast<float>(M(i, 0));
        om.r[i].y = static_cast<float>(M(i, 1));
        om.r[i].z = static_cast<float>(M(i, 2));
    }
    return om;
}

/// Score how well MOSFLM AMAT indexes reciprocal xyz (x = A · hkl).
IndexFitStats score_amat_vs_xyz(
    const Eigen::Matrix3d& amat,
    const std::vector<float>& x,
    const std::vector<float>& y,
    const std::vector<float>& z,
    float hkldev
) {
    IndexFitStats stats;
    const std::size_t n = x.size();
    stats.n_total = n;
    if (n == 0) {
        return stats;
    }
    Eigen::FullPivLU<Eigen::Matrix3d> lu(amat);
    if (!lu.isInvertible()) {
        return stats;
    }
    const Eigen::Matrix3d Ainv = lu.inverse();
    const double tol = static_cast<double>(hkldev);
    double sum_sq = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const Eigen::Vector3d r(
            static_cast<double>(x[i]),
            static_cast<double>(y[i]),
            static_cast<double>(z[i])
        );
        const Eigen::Vector3d hkl = Ainv * r;
        double max_dev = 0.0;
        double sq = 0.0;
        for (int j = 0; j < 3; ++j) {
            const double nearest = std::round(hkl[j]);
            const double d = hkl[j] - nearest;
            max_dev = std::max(max_dev, std::fabs(d));
            sq += d * d;
        }
        if (max_dev <= tol) {
            ++stats.n_indexed;
            sum_sq += sq;
        }
    }
    stats.indexed_frac =
        static_cast<double>(stats.n_indexed) / static_cast<double>(n);
    if (stats.n_indexed > 0) {
        stats.rmsd_hkl =
            std::sqrt(sum_sq / (3.0 * static_cast<double>(stats.n_indexed)));
    }
    return stats;
}

bool fit_passes(
    const IndexFitStats& fit,
    float min_frac,
    float max_rmsd
) {
    if (fit.n_total == 0) {
        return false;
    }
    if (fit.indexed_frac + 1e-12 < static_cast<double>(min_frac)) {
        return false;
    }
    if (fit.n_indexed > 0 &&
        fit.rmsd_hkl > static_cast<double>(max_rmsd) + 1e-12) {
        return false;
    }
    return true;
}

/// Prefer primitive monoclinic over C-centered when fits are comparable.
int monoclinic_prefer_key(std::string_view system_id) {
    if (system_id == "mP") {
        return 2;
    }
    if (system_id == "mC" || system_id == "mI") {
        return 0;
    }
    return 1;
}

}  // namespace

AutoIndexResult AutoIndexer::index(
    const SpotList& spots,
    const AutoIndexGeometry& geometry,
    const AutoIndexSettings& settings
) const {
    AutoIndexResult result;
    if (spots.spots.empty()) {
        return result;
    }

    const double beam_x =
        (geometry.beam_x_mm > 0.0) ? geometry.beam_x_mm : spots.header.beam_x_mm;
    const double beam_y =
        (geometry.beam_y_mm > 0.0) ? geometry.beam_y_mm : spots.header.beam_y_mm;

    std::vector<std::size_t> keep;
    keep.reserve(spots.spots.size());
    for (std::size_t i = 0; i < spots.spots.size(); ++i) {
        const auto& s = spots.spots[i];
        const double sig = (s.sigma > 1e-6) ? s.sigma : 1.0;
        if (s.intensity / sig >= static_cast<double>(settings.min_isigi)) {
            keep.push_back(i);
        }
    }
    const int n = static_cast<int>(keep.size());
    if (n > N_RS) {
        throw std::runtime_error("AutoIndexer: too many spots for DPS N_RS");
    }
    if (n < 20) {
        result.ierr = -2;
        result.n_spots_used = static_cast<std::size_t>(n);
        return result;
    }

    // MOSFLM DPSCHOOSE: FR_MAT (OMEGA + INVERTX) → triangular xyz → ROTMATY.
    // .spt INVFLAG≠0 ⇒ FINVERTX=-1 (ADSC). ROTMATY takes φ in degrees (DTOR).
    // Spot-line field 4 is φ_end (insulin: 0.5°). At TWOTHETA=MU=0 the γ/ν
    // path equals this triangular form.
    constexpr double kDeg2Rad = 3.14159265358979323846 / 180.0;
    const double dist = geometry.distance_mm;
    const double inv_l = 1.0 / geometry.wavelength_A;
    const double omega = spots.header.omega_deg * kDeg2Rad;
    const double finvertx = (spots.header.invflag != 0) ? -1.0 : 1.0;
    const double fr11 = std::cos(omega) * finvertx;
    const double fr12 = std::sin(omega);
    const double fr21 = -std::sin(omega) * finvertx;
    const double fr22 = std::cos(omega);
    const int n_xyz = n;
    result.n_spots_used = static_cast<std::size_t>(n_xyz);

    std::vector<float> x(static_cast<std::size_t>(n_xyz));
    std::vector<float> y(static_cast<std::size_t>(n_xyz));
    std::vector<float> z(static_cast<std::size_t>(n_xyz));
    std::vector<int> usespot(static_cast<std::size_t>(n_xyz), 1);
    for (int i = 0; i < n_xyz; ++i) {
        const auto& s = spots.spots[keep[static_cast<std::size_t>(i)]];
        const double xt = s.x_mm - beam_x;
        const double yt = s.y_mm - beam_y;
        const double xs = xt * fr11 + yt * fr12;
        const double ys = xt * fr21 + yt * fr22;
        const double rlen = std::sqrt(xs * xs + ys * ys + dist * dist);
        const double f = inv_l / rlen;
        const double qx = f * xs;
        const double qy = f * ys;
        const double qz = f * dist - inv_l;
        const double phi_rad = s.phi_end_deg * kDeg2Rad;
        const double c = std::cos(phi_rad);
        const double sn = std::sin(phi_rad);
        // ROTMATY (TRANSP∘RTMATS about Y): x'=c x − s z, z'=s x + c z
        x[static_cast<std::size_t>(i)] = static_cast<float>(c * qx - sn * qz);
        y[static_cast<std::size_t>(i)] = static_cast<float>(qy);
        z[static_cast<std::size_t>(i)] = static_cast<float>(sn * qx + c * qz);
    }

    int numvect = settings.nvect;
    int numvectref = settings.nvref;
    int maxcyc = settings.maxcyc;
    int inicyc = settings.inicyc;
    float convrg = settings.convrg;
    int numormat = settings.numormat;
    float stepscal = settings.stepscal;
    int numrlp = settings.numrlp;
    int numsptdps = settings.numsptdps;
    int npkfftmin = settings.npkfftmin;
    int nmultfft = settings.nmultfft;
    int idumpth = 0;
    float theta1d = settings.theta1d;
    float theta2d = settings.theta2d;
    float sinthmin = settings.sinthmin;
    float volmin = settings.volmin;
    float hkldev = settings.hkldev;
    float hklmin = settings.hklmin;
    int dumpfft[3] = {0, 0, 0};
    int dumpvec = 0;
    int dumpref = 0;
    int dumpgen = 0;
    int dumpch = 0;
    int number_xyz = n_xyz;
    int nlattice = 1;
    float d_max = static_cast<float>(geometry.dmax_A * settings.dmaxscal);
    float d_min = 0.0f;  // UMINCELL
    // MOSFLM DPSIND passes literal 1 (INDEX_AUTOMATIC is 0; unused in do_index).
    int index_mode = 1;
    float a9[9] = {};
    float wave = static_cast<float>(geometry.wavelength_A);
    int ierr = 0;
    float tarcell[6] = {0, 0, 0, 0, 0, 0};
    int centering = 0;

    index_fortran_(
        &numvect,
        &numvectref,
        &maxcyc,
        &inicyc,
        &convrg,
        &numormat,
        &stepscal,
        &numrlp,
        &numsptdps,
        &npkfftmin,
        &nmultfft,
        &idumpth,
        &theta1d,
        &theta2d,
        &sinthmin,
        &volmin,
        &hkldev,
        &hklmin,
        dumpfft,
        &dumpvec,
        &dumpref,
        &dumpgen,
        &dumpch,
        &number_xyz,
        &nlattice,
        x.data(),
        y.data(),
        z.data(),
        usespot.data(),
        &d_max,
        &d_min,
        &index_mode,
        a9,
        &wave,
        &ierr,
        tarcell,
        &centering
    );

    result.ierr = ierr;
    if (ierr != 0) {
        return result;
    }

    result.amat << a9[0], a9[1], a9[2], a9[3], a9[4], a9[5], a9[6], a9[7],
        a9[8];

    // Recover DPS OM from AMAT packing for cell + Bravais.
    struct matrix_3x3 om{};
    const double inv_w = 1.0 / geometry.wavelength_A;
    om.r[0].z = static_cast<float>(a9[0] * inv_w);
    om.r[0].x = static_cast<float>(a9[1] * inv_w);
    om.r[0].y = static_cast<float>(a9[2] * inv_w);
    om.r[1].z = static_cast<float>(a9[3] * inv_w);
    om.r[1].x = static_cast<float>(a9[4] * inv_w);
    om.r[1].y = static_cast<float>(a9[5] * inv_w);
    om.r[2].z = static_cast<float>(a9[6] * inv_w);
    om.r[2].x = static_cast<float>(a9[7] * inv_w);
    om.r[2].y = static_cast<float>(a9[8] * inv_w);

    struct cell_struct cell{};
    dircell_from_om(om, &cell);
    result.cell = cell_from_dps(cell);

    Eigen::Matrix3d om_e;
    for (int i = 0; i < 3; ++i) {
        om_e(i, 0) = om.r[i].x;
        om_e(i, 1) = om.r[i].y;
        om_e(i, 2) = om.r[i].z;
    }
    result.bravais = BravaisClassifier{}.classify(om_e, 1);

    const std::size_t sugg = select_suggested_bravais(
        result.bravais,
        settings.max_penalty,
        settings.prefer_character,
        settings.prefer_system
    );
    const bool user_prefer =
        settings.prefer_character > 0 || !settings.prefer_system.empty();
    const int lastgood = find_lastgood_rank(result.bravais);

    struct TrialScore {
        std::size_t idx{0};
        IndexFitStats fit{};
        Eigen::Matrix3d amat{Eigen::Matrix3d::Identity()};
        UnitCell cell{};
    };

    auto apply_candidate = [&](std::size_t idx) -> TrialScore {
        TrialScore ts;
        ts.idx = idx;
        if (idx >= result.bravais.size()) {
            return ts;
        }
        const auto& c = result.bravais[idx];
        const struct matrix_3x3 om_pick = eigen_to_dps(c.orientation);
        ts.amat = amat_from_dps_om(om_pick, geometry.wavelength_A);
        struct cell_struct cell2{};
        dircell_from_om(om_pick, &cell2);
        ts.cell = cell_from_dps(cell2);
        ts.fit = score_amat_vs_xyz(ts.amat, x, y, z, settings.hkldev);
        return ts;
    };

    std::size_t chosen = static_cast<std::size_t>(-1);
    TrialScore chosen_ts;

    if (user_prefer && sugg != static_cast<std::size_t>(-1) &&
        sugg < result.bravais.size()) {
        // Explicit SYSTEM / CHARACTER: keep SUGGLATT prefer pick (insulin path).
        chosen_ts = apply_candidate(sugg);
        chosen = sugg;
    } else {
        // Unconstrained: BESTSOLAR-lite among LASTGOOD (+ aP fallback).
        std::vector<std::size_t> trial;
        trial.reserve(result.bravais.size());
        auto push_unique = [&](std::size_t idx) {
            if (idx >= result.bravais.size()) {
                return;
            }
            if (std::find(trial.begin(), trial.end(), idx) == trial.end()) {
                trial.push_back(idx);
            }
        };
        push_unique(sugg);
        std::vector<std::size_t> by_pen(result.bravais.size());
        for (std::size_t i = 0; i < by_pen.size(); ++i) {
            by_pen[i] = i;
        }
        std::stable_sort(
            by_pen.begin(), by_pen.end(),
            [&](std::size_t a, std::size_t b) {
                if (result.bravais[a].criterion !=
                    result.bravais[b].criterion) {
                    return result.bravais[a].criterion <
                           result.bravais[b].criterion;
                }
                return result.bravais[a].character <
                       result.bravais[b].character;
            }
        );
        const int lim =
            std::max(1, std::min(lastgood, static_cast<int>(by_pen.size())));
        for (int rank = 1; rank <= lim; ++rank) {
            push_unique(by_pen[static_cast<std::size_t>(rank - 1)]);
        }
        for (std::size_t i = 0; i < result.bravais.size(); ++i) {
            if (result.bravais[i].system_id == "aP") {
                push_unique(i);
            }
        }

        std::vector<TrialScore> scored;
        scored.reserve(trial.size());
        for (const std::size_t idx : trial) {
            scored.push_back(apply_candidate(idx));
        }

        // Prefer mP over mC when indexed fractions are within 5%.
        std::stable_sort(
            scored.begin(), scored.end(),
            [&](const TrialScore& a, const TrialScore& b) {
                const auto& ca = result.bravais[a.idx];
                const auto& cb = result.bravais[b.idx];
                const double fa = a.fit.indexed_frac;
                const double fb = b.fit.indexed_frac;
                if (std::fabs(fa - fb) > 0.05) {
                    return fa > fb;
                }
                const int pa = monoclinic_prefer_key(ca.system_id);
                const int pb = monoclinic_prefer_key(cb.system_id);
                if (pa != pb) {
                    return pa > pb;
                }
                if (a.fit.rmsd_hkl != b.fit.rmsd_hkl) {
                    return a.fit.rmsd_hkl < b.fit.rmsd_hkl;
                }
                return ca.criterion < cb.criterion;
            }
        );

        for (const auto& ts : scored) {
            if (fit_passes(
                    ts.fit, settings.min_indexed_frac, settings.max_rmsd_hkl
                )) {
                chosen = ts.idx;
                chosen_ts = ts;
                break;
            }
        }
        // No gate passer: keep SUGGLATT unless an mP clearly out-indexes mC.
        if (chosen == static_cast<std::size_t>(-1) && !scored.empty()) {
            TrialScore sugg_ts;
            bool have_sugg = false;
            for (const auto& ts : scored) {
                if (ts.idx == sugg) {
                    sugg_ts = ts;
                    have_sugg = true;
                    break;
                }
            }
            const auto& best = scored.front();
            const bool sugg_is_mc =
                have_sugg &&
                (result.bravais[sugg].system_id == "mC" ||
                 result.bravais[sugg].system_id == "mI");
            const bool best_is_mp =
                result.bravais[best.idx].system_id == "mP";
            if (sugg_is_mc && best_is_mp &&
                best.fit.indexed_frac + 1e-12 >=
                    sugg_ts.fit.indexed_frac - 0.02) {
                chosen = best.idx;
                chosen_ts = best;
            } else if (have_sugg) {
                chosen = sugg;
                chosen_ts = sugg_ts;
            } else {
                chosen = best.idx;
                chosen_ts = best;
            }
        }
    }

    if (chosen != static_cast<std::size_t>(-1) &&
        chosen < result.bravais.size()) {
        const auto& c = result.bravais[chosen];
        result.amat = chosen_ts.amat;
        result.cell = chosen_ts.cell;
        result.fit = chosen_ts.fit;
        result.selected_character = c.character;
        result.selected_system = c.system_id;
        result.selected_symmetry_number =
            symmetry_number_for_lattice(c.system_id);
        std::cout << std::fixed;
        std::cout << "  INDEXFIT: system=" << c.system_id
                  << " char=" << c.character << " indexed="
                  << static_cast<int>(
                         std::lround(100.0 * result.fit.indexed_frac)
                     )
                  << "% rmsd_hkl=" << std::setprecision(3)
                  << result.fit.rmsd_hkl << " (n=" << result.fit.n_indexed
                  << "/" << result.fit.n_total << ")\n";
    }

    result.orientation =
        OrientationMatrix::from_mosflm_orientation_matrix(result.amat);
    result.ok = result.cell.is_valid();
    return result;
}

}  // namespace mosflm::indexing
