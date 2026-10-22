#include "mosflm/refinement/PositionalRefinement.hpp"

#include "mosflm/geometry/MetricTensorProjector.hpp"
#include "mosflm/refinement/PostRefinement.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace mosflm::refinement {
namespace {

constexpr double kDeg = 3.14159265358979323846 / 180.0;

Eigen::Matrix3d rotation_x_deg(double a) {
    const double c = std::cos(a * kDeg);
    const double s = std::sin(a * kDeg);
    Eigen::Matrix3d R;
    R << 1.0, 0.0, 0.0, 0.0, c, -s, 0.0, s, c;
    return R;
}

Eigen::Matrix3d rotation_y_deg(double a) {
    const double c = std::cos(a * kDeg);
    const double s = std::sin(a * kDeg);
    Eigen::Matrix3d R;
    R << c, 0.0, s, 0.0, 1.0, 0.0, -s, 0.0, c;
    return R;
}

Eigen::Matrix3d rotation_z_deg(double a) {
    const double c = std::cos(a * kDeg);
    const double s = std::sin(a * kDeg);
    Eigen::Matrix3d R;
    R << c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0;
    return R;
}

/// Reciprocal orthogonalization ``A^{-T}`` for the cell of metric ``G``.
Eigen::Matrix3d bstar_of(const Eigen::Matrix3d& G) {
    const UnitCell cell = cell_from_metric_tensor(G);
    return cell.orthogonalization_matrix().inverse().transpose();
}

/// Closest orthogonal matrix to ``M`` (polar decomposition). Handedness is
/// kept: the MOSFLM → predictor frame change can make ``UB`` improper.
Eigen::Matrix3d nearest_orthogonal(const Eigen::Matrix3d& M) {
    const Eigen::JacobiSVD<Eigen::Matrix3d> svd(
        M, Eigen::ComputeFullU | Eigen::ComputeFullV
    );
    return svd.matrixU() * svd.matrixV().transpose();
}

/// Parameter layout: [rot_x, rot_y, rot_z, metric c_1..c_k, distance,
/// beam_x, beam_y, tilt_mrad, twist_mrad].
struct Model {
    Eigen::Matrix3d U0;
    const MetricTensorProjector* projector;
    int n_metric;

    [[nodiscard]] OrientationMatrix orientation(const Eigen::VectorXd& p) const {
        const Eigen::Matrix3d G =
            projector->metric_from_coefficients(p.segment(3, n_metric));
        const Eigen::Matrix3d R =
            rotation_z_deg(p(2)) * rotation_y_deg(p(1)) * rotation_x_deg(p(0));
        return OrientationMatrix(R * U0 * bstar_of(G));
    }
};

/// Residual vector (3 per observation: X, Y px and weighted φ; zero for
/// unused observations); throws on a degenerate metric.
Eigen::VectorXd residuals(
    const Model& model,
    const Eigen::VectorXd& p,
    double wavelength_A,
    DetectorGeometry det,
    const std::vector<PositionalObservation>& obs,
    const std::vector<char>& use,
    double phi_weight
) {
    const OrientationMatrix om = model.orientation(p);
    const int base = 3 + model.n_metric;
    det.distance_mm = p(base);
    det.use_rotmat2 = true;
    det.tilt_mrad = p(base + 3);
    det.twist_mrad = p(base + 4);
    const double bx = p(base + 1);
    const double by = p(base + 2);
    Eigen::VectorXd r = Eigen::VectorXd::Zero(3 * static_cast<Eigen::Index>(obs.size()));
    for (std::size_t i = 0; i < obs.size(); ++i) {
        if (!use[i]) {
            continue;
        }
        const auto& o = obs[i];
        const bool has_phi = std::isfinite(o.phi_obs_deg);
        const double phi = nearest_ewald_phi_deg(
            om.miller_to_lab(o.h, o.k, o.l), wavelength_A,
            has_phi ? o.phi_obs_deg : o.phi_deg
        );
        double x = 0.0;
        double y = 0.0;
        const double w = std::sqrt(std::max(o.weight, 0.0));
        if (!std::isfinite(phi) ||
            !predict_xdet_ydet_px(om, o.h, o.k, o.l, phi, wavelength_A, det, bx, by, x, y)) {
            // Lost from the detector: large penalty keeps the step honest.
            r(3 * i) = 50.0 * w;
            r(3 * i + 1) = 50.0 * w;
            continue;
        }
        r(3 * i) = w * (x - o.xdet_obs_px);
        r(3 * i + 1) = w * (y - o.ydet_obs_px);
        if (has_phi) {
            r(3 * i + 2) = std::isfinite(o.phi_sigma_deg) && o.phi_sigma_deg > 0.0
                               ? w * (phi - o.phi_obs_deg) / o.phi_sigma_deg
                               : w * phi_weight * (phi - o.phi_obs_deg);
        }
    }
    return r;
}

double rms_of(const Eigen::VectorXd& r, const std::vector<char>& use) {
    const auto n = std::count(use.begin(), use.end(), char{1});
    return n > 0 ? std::sqrt(r.squaredNorm() / static_cast<double>(n)) : 0.0;
}

}  // namespace

PositionalRefinementResult PositionalRefiner::refine(
    const OrientationMatrix& orientation0,
    const gemmi::SpaceGroup* spacegroup,
    double wavelength_A,
    const DetectorGeometry& detector,
    double beam_x_mm,
    double beam_y_mm,
    const std::vector<PositionalObservation>& observations,
    const PositionalRefinementSettings& settings
) const {
    PositionalRefinementResult out;
    out.orientation = orientation0;
    out.distance_mm = detector.distance_mm;
    out.beam_x_mm = beam_x_mm;
    out.beam_y_mm = beam_y_mm;
    out.tilt_mrad = detector.tilt_mrad;
    out.twist_mrad = detector.twist_mrad;
    out.n_input = observations.size();
    if (observations.size() < settings.min_observations || !(wavelength_A > 0.0)) {
        out.reason = "too few observations";
        return out;
    }

    const MetricTensorProjector projector(
        spacegroup != nullptr ? spacegroup : gemmi::find_spacegroup_by_number(1)
    );
    // Direct metric from UB: G = (UBᵀ UB)⁻¹, projected onto the allowed set.
    const Eigen::Matrix3d ub0 = orientation0.ub();
    const Eigen::Matrix3d G_obs = (ub0.transpose() * ub0).inverse();
    const Eigen::VectorXd c0 = projector.coefficients(G_obs);
    const Eigen::Matrix3d G0 = projector.metric_from_coefficients(c0);
    Model model{
        nearest_orthogonal(ub0 * bstar_of(G0).inverse()), &projector,
        static_cast<int>(c0.size())
    };

    const int base = 3 + model.n_metric;
    const int np = base + 5;
    Eigen::VectorXd p = Eigen::VectorXd::Zero(np);
    p.segment(3, model.n_metric) = c0;
    const double phi_w = settings.phi_weight_px_per_deg;
    const bool any_phi = std::any_of(
        observations.begin(), observations.end(),
        [](const PositionalObservation& o) { return std::isfinite(o.phi_obs_deg); }
    );
    p(base) = detector.distance_mm;
    p(base + 1) = beam_x_mm;
    p(base + 2) = beam_y_mm;
    p(base + 3) = detector.tilt_mrad;
    p(base + 4) = detector.twist_mrad;
    const Eigen::VectorXd p_start = p;

    // Finite-difference steps and active mask.
    Eigen::VectorXd step(np);
    std::vector<char> active(static_cast<std::size_t>(np), 0);
    step(0) = step(1) = step(2) = 0.01;
    active[0] = settings.refine_orientation && any_phi;  // spindle: needs φ
    active[1] = active[2] = settings.refine_orientation;
    for (int k = 0; k < model.n_metric; ++k) {
        step(3 + k) = 1e-5 * std::max(1.0, std::fabs(c0(k)));
        active[static_cast<std::size_t>(3 + k)] = settings.refine_cell;
    }
    step(base) = 0.01;
    active[static_cast<std::size_t>(base)] = settings.refine_distance;
    step(base + 1) = step(base + 2) = 0.002;
    active[static_cast<std::size_t>(base + 1)] = settings.refine_beam;
    active[static_cast<std::size_t>(base + 2)] = settings.refine_beam;
    step(base + 3) = step(base + 4) = 0.1;
    active[static_cast<std::size_t>(base + 3)] = settings.refine_tilt_twist;
    active[static_cast<std::size_t>(base + 4)] = settings.refine_tilt_twist;

    std::vector<char> use(observations.size(), 1);
    const std::vector<char> all(observations.size(), 1);
    try {
        Eigen::VectorXd r = residuals(model, p, wavelength_A, detector, observations, use, phi_w);
        out.rms_initial_px = rms_of(r, use);
        double lambda = 1e-3;
        for (int iter = 0; iter < settings.max_iterations; ++iter) {
            // Outlier rejection: judge every spot against the current model
            // (masked residuals are zero, so recompute them unmasked).
            const double cut = std::max(
                settings.outlier_floor_px, settings.outlier_nsigma * rms_of(r, use)
            );
            const Eigen::VectorXd r_all =
                residuals(model, p, wavelength_A, detector, observations, all, phi_w);
            for (std::size_t i = 0; i < observations.size(); ++i) {
                const double w = std::sqrt(std::max(observations[i].weight, 1e-12));
                use[i] = r_all.segment(3 * static_cast<Eigen::Index>(i), 3).norm() / w <= cut ? 1 : 0;
            }
            if (static_cast<std::size_t>(std::count(use.begin(), use.end(), char{1})) <
                settings.min_observations) {
                out.reason = "too few observations after outlier rejection";
                return out;
            }
            r = residuals(model, p, wavelength_A, detector, observations, use, phi_w);
            const double cost = r.squaredNorm();

            Eigen::MatrixXd J = Eigen::MatrixXd::Zero(r.size(), np);
            for (int k = 0; k < np; ++k) {
                if (!active[static_cast<std::size_t>(k)]) {
                    continue;
                }
                Eigen::VectorXd q = p;
                q(k) += step(k);
                J.col(k) = (residuals(model, q, wavelength_A, detector, observations, use, phi_w) - r) /
                           step(k);
            }
            const Eigen::MatrixXd JtJ = J.transpose() * J;
            const Eigen::VectorXd g = J.transpose() * r;
            bool stepped = false;
            for (int tries = 0; tries < 8 && !stepped; ++tries) {
                Eigen::MatrixXd A = JtJ;
                for (int k = 0; k < np; ++k) {
                    A(k, k) = active[static_cast<std::size_t>(k)]
                                  ? JtJ(k, k) * (1.0 + lambda) + 1e-12
                                  : 1.0;
                }
                const Eigen::VectorXd dp = A.ldlt().solve(-g);
                const Eigen::VectorXd q = p + dp;
                const Eigen::VectorXd rq =
                    residuals(model, q, wavelength_A, detector, observations, use, phi_w);
                if (rq.squaredNorm() < cost) {
                    p = q;
                    r = rq;
                    lambda = std::max(lambda * 0.3, 1e-7);
                    stepped = true;
                } else {
                    lambda *= 10.0;
                }
            }
            if (!stepped || std::fabs(cost - r.squaredNorm()) < 1e-6 * cost) {
                break;
            }
        }
        out.rms_px = rms_of(r, use);
        // Compare like with like: starting residual on the final subset.
        out.rms_initial_px = rms_of(
            residuals(model, p_start, wavelength_A, detector, observations, use, phi_w), use
        );
        out.n_used = static_cast<std::size_t>(std::count(use.begin(), use.end(), char{1}));

        // Safety caps on the total change.
        const Eigen::VectorXd dpt = p - p_start;
        const double metric_frac =
            (projector.metric_from_coefficients(p.segment(3, model.n_metric)) - G0).norm() /
            G0.norm();
        if (std::fabs(dpt(0)) > settings.max_rotation_deg ||
            std::fabs(dpt(1)) > settings.max_rotation_deg ||
            std::fabs(dpt(2)) > settings.max_rotation_deg ||
            std::fabs(dpt(base)) > settings.max_distance_frac * detector.distance_mm ||
            std::hypot(dpt(base + 1), dpt(base + 2)) > settings.max_beam_shift_mm ||
            metric_frac > settings.max_metric_frac ||
            std::fabs(p(base + 3)) > settings.max_tilt_twist_mrad ||
            std::fabs(p(base + 4)) > settings.max_tilt_twist_mrad) {
            out.reason = "change exceeds caps: rot " + std::to_string(dpt(0)) + "," +
                         std::to_string(dpt(1)) + "," + std::to_string(dpt(2)) +
                         " deg, distance " +
                         std::to_string(dpt(base)) + " mm, beam " +
                         std::to_string(std::hypot(dpt(base + 1), dpt(base + 2))) +
                         " mm, metric " + std::to_string(metric_frac) +
                         ", tilt/twist " + std::to_string(p(base + 3)) + "/" +
                         std::to_string(p(base + 4)) + " mrad" +
                         ", rms " + std::to_string(rms_of(r, use)) + " px";
            return out;
        }
        out.orientation = model.orientation(p);
        out.cell = cell_from_metric_tensor(
            projector.metric_from_coefficients(p.segment(3, model.n_metric))
        );
        out.distance_mm = p(base);
        out.beam_x_mm = p(base + 1);
        out.beam_y_mm = p(base + 2);
        out.tilt_mrad = p(base + 3);
        out.twist_mrad = p(base + 4);
        out.ok = true;
        out.improved = out.rms_px < out.rms_initial_px;
    } catch (const std::domain_error& ex) {
        out.ok = false;
        out.reason = ex.what();
    }
    return out;
}

}  // namespace mosflm::refinement
