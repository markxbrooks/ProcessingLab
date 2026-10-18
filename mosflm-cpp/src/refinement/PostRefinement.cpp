#include "mosflm/refinement/PostRefinement.hpp"

#include "mosflm/geometry/Rotmat.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace mosflm::refinement {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kRad = kPi / 180.0;

Eigen::Matrix3d rotation_xyz_deg(double dx, double dy, double dz) {
    return rotmat_from_missetting(dx, dy, dz, true);
}

UnitCell scale_cell(const UnitCell& c0, double s) {
    UnitCell c = c0;
    c.a *= s;
    c.b *= s;
    c.c *= s;
    return c;
}

struct ResidualStats {
    double rms{0.0};
    std::size_t n{0};
    std::vector<double> abs_resid;
};

ResidualStats evaluate_residuals(
    const OrientationMatrix& orient,
    double wavelength_A,
    const std::vector<PostRefObservation>& obs,
    double outlier_phi_deg
) {
    ResidualStats st;
    double wsum = 0.0;
    double wres2 = 0.0;
    for (const auto& o : obs) {
        const Eigen::Vector3d q0 = orient.miller_to_lab(o.h, o.k, o.l);
        const double phi_pred =
            nearest_ewald_phi_deg(q0, wavelength_A, o.phi_obs_deg);
        if (!std::isfinite(phi_pred)) {
            continue;
        }
        const double d = phi_pred - o.phi_obs_deg;
        if (std::fabs(d) > outlier_phi_deg) {
            continue;
        }
        const double w = (o.weight > 0.0) ? o.weight : 1.0;
        wsum += w;
        wres2 += w * d * d;
        st.abs_resid.push_back(std::fabs(d));
        ++st.n;
    }
    st.rms = (wsum > 0.0) ? std::sqrt(wres2 / wsum) : 0.0;
    return st;
}

}  // namespace

std::vector<double> ewald_phi_solutions_deg(
    const Eigen::Vector3d& q0,
    double wavelength_A
) {
    std::vector<double> out;
    if (!(wavelength_A > 0.0)) {
        return out;
    }
    const double qx = q0.x();
    const double qy = q0.y();
    const double qz = q0.z();
    (void)qx;
    const double q2 = q0.squaredNorm();
    // qy sinφ + qz cosφ = −(λ/2) |q|²
    const double rhs = -0.5 * wavelength_A * q2;
    const double amp = std::hypot(qy, qz);
    if (amp < 1e-15) {
        return out;
    }
    if (std::fabs(rhs) > amp + 1e-12) {
        return out;
    }
    const double cos_alpha = qz / amp;  // amp * cos(φ − α) = rhs with
    const double sin_alpha = qy / amp;  // α = atan2(qy, qz)
    const double alpha = std::atan2(sin_alpha, cos_alpha);
    const double arg = std::clamp(rhs / amp, -1.0, 1.0);
    const double delta = std::acos(arg);
    const double phi1 = (alpha + delta) / kRad;
    const double phi2 = (alpha - delta) / kRad;
    out.push_back(phi1);
    if (std::fabs(delta) > 1e-9) {
        out.push_back(phi2);
    }
    // Normalize to (−180, 180]
    for (double& p : out) {
        while (p > 180.0) {
            p -= 360.0;
        }
        while (p <= -180.0) {
            p += 360.0;
        }
    }
    return out;
}

double nearest_ewald_phi_deg(
    const Eigen::Vector3d& q0,
    double wavelength_A,
    double target_deg
) {
    const auto sols = ewald_phi_solutions_deg(q0, wavelength_A);
    if (sols.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double best = sols[0];
    double best_d = std::fabs(sols[0] - target_deg);
    for (double sol : sols) {
        for (double wrap : {0.0, 360.0, -360.0}) {
            const double cand = sol + wrap;
            const double d = std::fabs(cand - target_deg);
            if (d < best_d) {
                best_d = d;
                best = cand;
            }
        }
    }
    return best;
}

OrientationMatrix PostRefiner::apply_parameters(
    const UnitCell& cell0,
    const OrientationMatrix& orientation0,
    const PostRefParameters& p
) {
    const UnitCell cell = scale_cell(cell0, p.cell_scale);
    // Extract rotation from UB0 ≈ R0 * A*(cell0); rebuild with scaled cell.
    // Approximate: UB_new = R_misset * UB0 * (A*_new / A*_old) via
    // UB_new = R_misset * UB0 * (cell0.a / cell.a) for isotropic scale of
    // reciprocal (A* scales as 1/a).
    const double recip_scale =
        (p.cell_scale > 1e-12) ? (1.0 / p.cell_scale) : 1.0;
    const Eigen::Matrix3d R =
        rotation_xyz_deg(p.dphix_deg, p.dphiy_deg, p.dphiz_deg);
    return OrientationMatrix(R * (orientation0.ub() * recip_scale));
}

PostRefinementResult PostRefiner::refine(
    const crystal::CrystalMetadata& crystal,
    const OrientationMatrix& orientation0,
    double wavelength_A,
    double mosaic0_deg,
    const std::vector<PostRefObservation>& observations,
    const PostRefinementSettings& settings
) const {
    PostRefinementResult result;
    result.n_input = observations.size();
    result.parameters.mosaic_spread_deg = mosaic0_deg;
    result.orientation = orientation0;
    result.cell = crystal.to_mosflm_cell();

    if (observations.size() < 8 || !(wavelength_A > 0.0)) {
        return result;
    }

    const UnitCell cell0 = result.cell;
    PostRefParameters p;
    p.mosaic_spread_deg = mosaic0_deg;

    auto orient_of = [&](const PostRefParameters& pp) {
        return apply_parameters(cell0, orientation0, pp);
    };

    ResidualStats st0 = evaluate_residuals(
        orientation0, wavelength_A, observations, settings.outlier_phi_deg
    );
    result.rms_phi_initial_deg = st0.rms;

    ResidualStats best_st = st0;
    PostRefParameters best_p = p;

    for (int iter = 0; iter < settings.max_iterations; ++iter) {
        // Coordinate descent with finite differences on active parameters.
        struct Step {
            double* value;
            double delta;
            bool active;
        };
        const Step steps[] = {
            {&p.dphix_deg, settings.fd_misset_deg, settings.refine_missetting},
            {&p.dphiy_deg, settings.fd_misset_deg, settings.refine_missetting},
            {&p.dphiz_deg, settings.fd_misset_deg, settings.refine_missetting},
            {&p.cell_scale, settings.fd_cell_scale, settings.refine_cell_scale},
        };

        bool improved = false;
        for (const Step& step : steps) {
            if (!step.active) {
                continue;
            }
            const double v0 = *step.value;
            double best_v = v0;
            double best_rms = best_st.rms;
            for (double s : {-1.0, 1.0, -2.0, 2.0}) {
                *step.value = v0 + s * step.delta;
                if (step.value == &p.cell_scale &&
                    (*step.value < 0.9 || *step.value > 1.1)) {
                    continue;
                }
                const ResidualStats st = evaluate_residuals(
                    orient_of(p), wavelength_A, observations,
                    settings.outlier_phi_deg
                );
                if (st.n >= 8 && st.rms < best_rms) {
                    best_rms = st.rms;
                    best_v = *step.value;
                    best_st = st;
                    improved = true;
                }
            }
            *step.value = best_v;
        }
        best_p = p;
        if (!improved) {
            break;
        }
    }

    // Mosaic from MAD of residuals (η ≈ 2.355 * σ_φ for Gaussian rocking
    // stand-in; clamp to a crystallographic range).
    if (settings.refine_mosaic && !best_st.abs_resid.empty()) {
        std::vector<double> r = best_st.abs_resid;
        std::sort(r.begin(), r.end());
        const double med = r[r.size() / 2];
        std::vector<double> dev;
        dev.reserve(r.size());
        for (double v : r) {
            dev.push_back(std::fabs(v - med));
        }
        std::sort(dev.begin(), dev.end());
        const double mad = dev[dev.size() / 2];
        const double sigma_phi = 1.4826 * mad;
        best_p.mosaic_spread_deg =
            std::clamp(2.355 * sigma_phi, 0.05, 2.0);
    }

    result.parameters = best_p;
    result.cell = scale_cell(cell0, best_p.cell_scale);
    result.orientation = orient_of(best_p);
    result.rms_phi_deg = best_st.rms;
    result.n_used = best_st.n;
    result.ok = best_st.n >= 8 && best_st.rms <= st0.rms + 1e-6;
    return result;
}

}  // namespace mosflm::refinement
