#include "mosflm/indexing/IdxrefRefine.hpp"

#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/geometry/LabFrame.hpp"
#include "mosflm/geometry/Rotmat.hpp"
#include "mosflm/refinement/PostRefinement.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace mosflm::indexing {
namespace {

Eigen::Matrix3d rotation_xyz_deg(double dx, double dy, double dz) {
    return rotmat_from_missetting(dx, dy, dz, true);
}

Eigen::Matrix3d astar_of(const UnitCell& cell) {
    const Eigen::Matrix3d A = cell.orthogonalization_matrix();
    return A.inverse().transpose();
}

struct Params {
    double dphix{0.0};
    double dphiy{0.0};
    double dphiz{0.0};
    UnitCell cell{};
    double beam_x{0.0};
    double beam_y{0.0};
    double distance_mm{0.0};
};

OrientationMatrix orient_of(
    const Eigen::Matrix3d& R0,
    const Params& p
) {
    const Eigen::Matrix3d Rmis =
        rotation_xyz_deg(p.dphix, p.dphiy, p.dphiz);
    return OrientationMatrix(Rmis * R0 * astar_of(p.cell));
}

struct ResidualStats {
    double rms_phi{0.0};
    double rms_pos{0.0};
    double cost{0.0};
    std::size_t n{0};
    std::size_t n_pos{0};
};

ResidualStats evaluate(
    const OrientationMatrix& orient,
    const Eigen::Matrix3d& amat_mosflm,
    const Params& p,
    double wavelength_A,
    double pixel_mm,
    double yscale,
    const std::vector<IdxrefObservation>& obs,
    const IdxrefSettings& settings
) {
    ResidualStats st;
    double wsum_phi = 0.0;
    double wres2_phi = 0.0;
    double wsum_pos = 0.0;
    double wres2_pos = 0.0;
    double wsum_hkl = 0.0;
    double wres2_hkl = 0.0;

    DetectorGeometry det;
    det.distance_mm = p.distance_mm;
    det.wavelength_A = wavelength_A;
    det.pixel_size_x_mm = pixel_mm;
    det.pixel_size_y_mm = pixel_mm;
    det.yscale = yscale;
    det.beam_x_mm = 0.0;
    det.beam_y_mm = 0.0;

    for (const auto& o : obs) {
        const double w = (o.weight > 0.0) ? o.weight : 1.0;
        bool counted = false;

        if (std::isfinite(o.rx) && std::isfinite(o.ry) && std::isfinite(o.rz)) {
            const Eigen::Vector3d r(o.rx, o.ry, o.rz);
            const Eigen::Vector3d pred =
                amat_mosflm * Eigen::Vector3d(o.h, o.k, o.l);
            const double dhkl = (r - pred).norm();
            wsum_hkl += w;
            wres2_hkl += w * dhkl * dhkl;
            counted = true;
        }

        const Eigen::Vector3d q0 = orient.miller_to_lab(o.h, o.k, o.l);
        const double phi_pred = refinement::nearest_ewald_phi_deg(
            q0, wavelength_A, o.phi_obs_deg
        );
        if (std::isfinite(phi_pred)) {
            const double dphi = phi_pred - o.phi_obs_deg;
            if (std::fabs(dphi) <= settings.outlier_phi_deg) {
                wsum_phi += w;
                wres2_phi += w * dphi * dphi;
                counted = true;

                double xdet = 0.0;
                double ydet = 0.0;
                if (refinement::predict_xdet_ydet_px(
                        orient, o.h, o.k, o.l, phi_pred, wavelength_A, det,
                        p.beam_x, p.beam_y, xdet, ydet
                    )) {
                    const double dpos = std::hypot(
                        xdet - o.xdet_obs_px, ydet - o.ydet_obs_px
                    );
                    if (dpos <= settings.outlier_pos_px) {
                        wsum_pos += w;
                        wres2_pos += w * dpos * dpos;
                        ++st.n_pos;
                    }
                }
            }
        }
        if (counted) {
            ++st.n;
        }
    }

    st.rms_phi = (wsum_phi > 0.0) ? std::sqrt(wres2_phi / wsum_phi) : 0.0;
    st.rms_pos = (wsum_pos > 0.0) ? std::sqrt(wres2_pos / wsum_pos) : 0.0;
    const double rms_hkl =
        (wsum_hkl > 0.0) ? std::sqrt(wres2_hkl / wsum_hkl) : 0.0;
    const double px_scale = std::max(settings.position_px_per_deg, 1.0);
    // Kabsch §2.8: primary XY + φ; HKL is a soft regularizer only.
    st.cost = st.rms_phi;
    if (st.n_pos > 0) {
        st.cost += st.rms_pos / px_scale;
    }
    st.cost += settings.hkl_cost_weight * rms_hkl;
    return st;
}

bool cell_ok(const UnitCell& c) {
    return c.is_valid() && c.a > 5.0 && c.b > 5.0 && c.c > 5.0 &&
           c.a < 2000.0 && c.b < 2000.0 && c.c < 2000.0 &&
           c.alpha_deg > 20.0 && c.beta_deg > 20.0 && c.gamma_deg > 20.0 &&
           c.alpha_deg < 160.0 && c.beta_deg < 160.0 && c.gamma_deg < 160.0;
}

double phi_span_deg(const std::vector<IdxrefObservation>& obs) {
    double lo = std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();
    for (const auto& o : obs) {
        if (!std::isfinite(o.phi_obs_deg)) {
            continue;
        }
        lo = std::min(lo, o.phi_obs_deg);
        hi = std::max(hi, o.phi_obs_deg);
    }
    if (!(hi > lo)) {
        return 0.0;
    }
    return hi - lo;
}

ResidualStats eval_params(
    const Eigen::Matrix3d& R0,
    const Params& p,
    double wavelength_A,
    double pixel_mm,
    double yscale,
    const std::vector<IdxrefObservation>& observations,
    const IdxrefSettings& settings
) {
    const OrientationMatrix o = orient_of(R0, p);
    const Eigen::Matrix3d am =
        cpp_ub_to_mosflm_amat(o.ub(), wavelength_A);
    return evaluate(
        o, am, p, wavelength_A, pixel_mm, yscale, observations, settings
    );
}

}  // namespace

std::vector<IdxrefObservation> build_idxref_observations(
    const SpotList& spots,
    const std::vector<std::size_t>& keep,
    const LocalIndexResult& local,
    const IdxrefGeometry& geometry,
    const std::vector<float>* x,
    const std::vector<float>* y,
    const std::vector<float>* z
) {
    std::vector<IdxrefObservation> out;
    out.reserve(keep.size());
    const double pix =
        (geometry.pixel_mm > 1e-9) ? geometry.pixel_mm : 0.075;
    for (std::size_t i = 0; i < keep.size(); ++i) {
        if (i >= local.spots.size() || !local.spots[i].indexed) {
            continue;
        }
        const auto& s = spots.spots[keep[i]];
        IdxrefObservation o;
        o.h = local.spots[i].h;
        o.k = local.spots[i].k;
        o.l = local.spots[i].l;
        o.phi_obs_deg = s.phi_end_deg;
        // .spt stores absolute image-frame mm (lower-left). Map to MOSFLM
        // XDET/YDET via the same helper as predict_xdet_ydet_px: undo the
        // GENDATA XG↔y / YG↔x swap so image (X,Y) → pixels match predictions.
        const double x_cpp = s.y_mm - geometry.beam_y_mm;
        const double y_cpp = s.x_mm - geometry.beam_x_mm;
        cpp_mm_to_mosflm_xdet_yscale(
            x_cpp,
            y_cpp,
            geometry.beam_x_mm,
            geometry.beam_y_mm,
            pix,
            pix,
            geometry.yscale,
            o.xdet_obs_px,
            o.ydet_obs_px
        );
        const double sig = (s.sigma > 1e-6) ? s.sigma : 1.0;
        o.weight = std::max(0.1, s.intensity / sig);
        if (x && y && z && i < x->size() && i < y->size() && i < z->size()) {
            o.rx = static_cast<double>((*x)[i]);
            o.ry = static_cast<double>((*y)[i]);
            o.rz = static_cast<double>((*z)[i]);
        }
        out.push_back(o);
    }
    return out;
}

Eigen::Vector3i select_index_origin(
    const UnitCell& cell0,
    const Eigen::Matrix3d& amat0_mosflm,
    std::vector<IdxrefObservation>& observations,
    const IdxrefGeometry& /*geometry*/,
    const IdxrefSettings& settings
) {
    Eigen::Vector3i best(0, 0, 0);
    if (observations.size() < 12 || !cell0.is_valid()) {
        return best;
    }

    // Kabsch §2.7: choose constant offset so observed RLPs and A·hkl
    // centroids are closest (HKL residual). Detector XY is unreliable here
    // before geometry LS and must not drive the origin pick.
    std::size_t n_xyz = 0;
    for (const auto& o : observations) {
        if (std::isfinite(o.rx) && std::isfinite(o.ry) && std::isfinite(o.rz)) {
            ++n_xyz;
        }
    }
    if (n_xyz < 8) {
        return best;
    }

    std::vector<Eigen::Vector3i> base;
    base.reserve(observations.size());
    for (const auto& o : observations) {
        base.emplace_back(o.h, o.k, o.l);
    }

    auto hkl_rms = [&](int dh, int dk, int dl) -> double {
        double wsum = 0.0;
        double wres2 = 0.0;
        for (std::size_t i = 0; i < observations.size(); ++i) {
            const auto& o = observations[i];
            if (!std::isfinite(o.rx) || !std::isfinite(o.ry) ||
                !std::isfinite(o.rz)) {
                continue;
            }
            const Eigen::Vector3d r(o.rx, o.ry, o.rz);
            const Eigen::Vector3d pred =
                amat0_mosflm * Eigen::Vector3d(
                    base[i].x() + dh, base[i].y() + dk, base[i].z() + dl
                );
            const double w = (o.weight > 0.0) ? o.weight : 1.0;
            const double d = (r - pred).norm();
            wsum += w;
            wres2 += w * d * d;
        }
        if (wsum < 1e-12) {
            return std::numeric_limits<double>::infinity();
        }
        return std::sqrt(wres2 / wsum);
    };

    const int n = std::max(0, settings.index_origin_max);
    double best_rms = hkl_rms(0, 0, 0);
    for (int dh = -n; dh <= n; ++dh) {
        for (int dk = -n; dk <= n; ++dk) {
            for (int dl = -n; dl <= n; ++dl) {
                if (dh == 0 && dk == 0 && dl == 0) {
                    continue;
                }
                const double rms = hkl_rms(dh, dk, dl);
                if (rms < best_rms - 1e-12) {
                    best_rms = rms;
                    best = Eigen::Vector3i(dh, dk, dl);
                }
            }
        }
    }

    for (std::size_t i = 0; i < observations.size(); ++i) {
        observations[i].h = base[i].x() + best.x();
        observations[i].k = base[i].y() + best.y();
        observations[i].l = base[i].z() + best.z();
    }
    return best;
}

IdxrefResult refine_idxref(
    const UnitCell& cell0,
    const Eigen::Matrix3d& amat0_mosflm,
    const std::vector<IdxrefObservation>& observations,
    const IdxrefGeometry& geometry,
    const IdxrefSettings& settings
) {
    IdxrefResult result;
    result.n_input = observations.size();
    result.cell = cell0;
    result.amat = amat0_mosflm;
    result.beam_x_mm = geometry.beam_x_mm;
    result.beam_y_mm = geometry.beam_y_mm;
    result.distance_mm = geometry.distance_mm;
    result.distance_mm_initial = geometry.distance_mm;
    result.orientation = OrientationMatrix::from_mosflm_orientation_matrix(
        amat0_mosflm, geometry.wavelength_A
    );

    if (observations.size() < 12 || !(geometry.wavelength_A > 0.0) ||
        !cell0.is_valid()) {
        return result;
    }

    const Eigen::Matrix3d Astar0 = astar_of(cell0);
    Eigen::FullPivLU<Eigen::Matrix3d> lu(Astar0);
    if (!lu.isInvertible()) {
        return result;
    }
    const Eigen::Matrix3d R0 =
        result.orientation.ub() * lu.inverse();

    Params p;
    p.cell = cell0;
    p.beam_x = geometry.beam_x_mm;
    p.beam_y = geometry.beam_y_mm;
    p.distance_mm = geometry.distance_mm;

    const bool refine_f =
        phi_span_deg(observations) >=
        settings.refine_distance_min_phi_span_deg;

    ResidualStats st0 = eval_params(
        R0, p, geometry.wavelength_A, geometry.pixel_mm, geometry.yscale,
        observations, settings
    );
    result.rms_phi_initial = st0.rms_phi;
    result.rms_pos_initial = st0.rms_pos;
    if (st0.n < 8) {
        return result;
    }

    ResidualStats best_st = st0;
    Params best_p = p;

    auto try_step = [&](double* value, double delta, double lo, double hi) {
        const double v0 = *value;
        double best_v = v0;
        double best_cost = best_st.cost;
        bool improved = false;
        for (double s : {-1.0, 1.0, -2.0, 2.0}) {
            *value = v0 + s * delta;
            if (*value < lo || *value > hi) {
                continue;
            }
            if ((value == &p.dphix || value == &p.dphiy ||
                 value == &p.dphiz) &&
                std::fabs(*value) > settings.max_misset_step_deg * 4.0) {
                continue;
            }
            if (!cell_ok(p.cell)) {
                continue;
            }
            const ResidualStats st = eval_params(
                R0, p, geometry.wavelength_A, geometry.pixel_mm,
                geometry.yscale, observations, settings
            );
            if (st.n >= 8 && st.cost < best_cost - 1e-9) {
                best_cost = st.cost;
                best_v = *value;
                best_st = st;
                improved = true;
            }
        }
        *value = best_v;
        return improved;
    };

    const double d0 = geometry.distance_mm;
    const double d_lo = d0 * (1.0 - settings.distance_frac_bound);
    const double d_hi = d0 * (1.0 + settings.distance_frac_bound);

    for (int iter = 0; iter < settings.max_iterations; ++iter) {
        bool improved = false;
        improved |= try_step(
            &p.dphix, settings.fd_misset_deg, -2.0, 2.0
        );
        improved |= try_step(
            &p.dphiy, settings.fd_misset_deg, -2.0, 2.0
        );
        improved |= try_step(
            &p.dphiz, settings.fd_misset_deg, -2.0, 2.0
        );
        improved |= try_step(
            &p.cell.a, settings.fd_cell_len, 0.85 * cell0.a, 1.15 * cell0.a
        );
        improved |= try_step(
            &p.cell.b, settings.fd_cell_len, 0.85 * cell0.b, 1.15 * cell0.b
        );
        improved |= try_step(
            &p.cell.c, settings.fd_cell_len, 0.85 * cell0.c, 1.15 * cell0.c
        );
        improved |= try_step(
            &p.cell.alpha_deg, settings.fd_cell_ang,
            cell0.alpha_deg - 8.0, cell0.alpha_deg + 8.0
        );
        improved |= try_step(
            &p.cell.beta_deg, settings.fd_cell_ang, cell0.beta_deg - 8.0,
            cell0.beta_deg + 8.0
        );
        improved |= try_step(
            &p.cell.gamma_deg, settings.fd_cell_ang,
            cell0.gamma_deg - 8.0, cell0.gamma_deg + 8.0
        );
        improved |= try_step(
            &p.beam_x, settings.fd_beam_mm, geometry.beam_x_mm - 2.0,
            geometry.beam_x_mm + 2.0
        );
        improved |= try_step(
            &p.beam_y, settings.fd_beam_mm, geometry.beam_y_mm - 2.0,
            geometry.beam_y_mm + 2.0
        );
        if (refine_f) {
            // Smaller relative steps to limit beam–distance coupling.
            improved |= try_step(
                &p.distance_mm, settings.fd_distance_mm, d_lo, d_hi
            );
        }
        best_p = p;
        if (!improved) {
            break;
        }
    }

    const OrientationMatrix orient_best = orient_of(R0, best_p);
    result.orientation = orient_best;
    result.cell = best_p.cell;
    result.beam_x_mm = best_p.beam_x;
    result.beam_y_mm = best_p.beam_y;
    result.distance_mm = best_p.distance_mm;
    result.amat = cpp_ub_to_mosflm_amat(
        orient_best.ub(), geometry.wavelength_A
    );
    result.rms_phi = best_st.rms_phi;
    result.rms_pos = best_st.rms_pos;
    result.n_used = best_st.n;
    result.ok = best_st.n >= 8;

    const bool pos_ok =
        best_st.rms_pos <= st0.rms_pos * 1.02 + 1e-9 ||
        st0.rms_pos <= 0.0;
    const bool phi_ok =
        best_st.rms_phi <= st0.rms_phi * 1.02 + 1e-9;
    result.accepted =
        result.ok && best_st.cost <= st0.cost + 1e-6 && (pos_ok || phi_ok);
    return result;
}

}  // namespace mosflm::indexing
