/**
 * PositionalRefiner: recover distance, beam, tetragonal cell and misset
 * from synthetic thaumatin-like spot centroids.
 */
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/refinement/PositionalRefinement.hpp"
#include "mosflm/refinement/PostRefinement.hpp"
#include "test_util.hpp"

#include <gemmi/symmetry.hpp>

#include <cmath>
#include <iostream>
#include <random>

namespace {

Eigen::Matrix3d rot_y(double deg) {
    const double a = deg * 3.14159265358979323846 / 180.0;
    Eigen::Matrix3d R;
    R << std::cos(a), 0, std::sin(a), 0, 1, 0, -std::sin(a), 0, std::cos(a);
    return R;
}

}  // namespace

int main() {
    using mosflm::OrientationMatrix;
    using mosflm::UnitCell;
    using namespace mosflm::refinement;

    const double lam = 0.9795;
    const auto* p4 = gemmi::find_spacegroup_by_number(75);
    const UnitCell cell_true{57.94, 57.94, 150.27, 90, 90, 90};
    const Eigen::Matrix3d U = Eigen::AngleAxisd(0.4, Eigen::Vector3d(1, 2, 3).normalized())
                                  .toRotationMatrix();
    const auto om_true = OrientationMatrix::from_cell_and_rotation(cell_true, U);

    mosflm::DetectorGeometry det_true;
    det_true.distance_mm = 199.78;
    det_true.wavelength_A = lam;
    det_true.pixel_size_x_mm = det_true.pixel_size_y_mm = 0.075;
    det_true.use_rotmat2 = true;
    det_true.tilt_mrad = 7.3;
    det_true.twist_mrad = -9.6;
    const double bx_true = 166.15, by_true = 172.97;

    std::mt19937 rng(5);
    std::uniform_int_distribution<int> hk(-30, 30), ll(-80, 80);
    std::normal_distribution<double> noise(0.0, 0.3);
    std::vector<PositionalObservation> obs;
    while (obs.size() < 600) {
        const int h = hk(rng), k = hk(rng), l = ll(rng);
        const Eigen::Vector3d q = om_true.miller_to_lab(h, k, l);
        if (q.norm() < 1e-9 || q.norm() > 1.0 / 1.6) {
            continue;
        }
        for (double phi : ewald_phi_solutions_deg(q, lam)) {
            if (phi < 0.0 || phi > 30.0) {
                continue;
            }
            double x = 0, y = 0;
            if (predict_xdet_ydet_px(om_true, h, k, l, phi, lam, det_true, bx_true, by_true, x, y) &&
                x > 0 && y > 0 && x < 4148 && y < 4362) {
                obs.push_back({h, k, l, phi + 0.05, x + noise(rng), y + noise(rng), 1.0});
            }
        }
    }
    // A few gross outliers.
    for (int i = 0; i < 20; ++i) {
        obs[static_cast<std::size_t>(i)].xdet_obs_px += 40.0;
    }

    // Start: CELLFIX-like cell (+0.27%), nominal distance, beam off, misset.
    const UnitCell cell0{58.09, 58.24, 150.33, 90.2, 90.1, 90.0};
    const auto om0 = OrientationMatrix::from_cell_and_rotation(cell0, rot_y(0.15) * U);
    mosflm::DetectorGeometry det0 = det_true;
    det0.distance_mm = 200.0;
    det0.tilt_mrad = 0.0;
    det0.twist_mrad = 0.0;

    const auto res = PositionalRefiner{}.refine(
        om0, p4, lam, det0, bx_true + 0.05, by_true - 0.04, obs
    );
    std::cout << "rms " << res.rms_initial_px << " -> " << res.rms_px << " n=" << res.n_used
              << " dist=" << res.distance_mm << " a=" << res.cell.a << " b=" << res.cell.b
              << " c=" << res.cell.c << " beam=" << res.beam_x_mm << "," << res.beam_y_mm << " tilt/twist=" << res.tilt_mrad
              << "/" << res.twist_mrad << '\n';
    CHECK(res.ok && res.improved);
    CHECK(res.rms_px < 0.5);
    CHECK(res.n_used >= 560u && res.n_used <= 590u);
    CHECK(std::fabs(res.distance_mm - 199.78) < 0.1);
    CHECK(std::fabs(res.cell.a - res.cell.b) < 1e-9);
    CHECK(std::fabs(res.cell.gamma_deg - 90.0) < 1e-9);
    CHECK(std::fabs(res.cell.a - 57.94) < 0.03);
    CHECK(std::fabs(res.cell.c - 150.27) < 0.1);
    CHECK(std::fabs(res.beam_x_mm - bx_true) < 0.01);
    CHECK(std::fabs(res.beam_y_mm - by_true) < 0.01);
    CHECK(std::fabs(res.tilt_mrad - 7.3) < 0.3);
    CHECK(std::fabs(res.twist_mrad + 9.6) < 0.3);

    // With observed φ, a misset about the spindle (+X) is recovered too:
    // without φ it is invisible in X/Y.
    {
        std::vector<PositionalObservation> obs_phi;
        std::normal_distribution<double> phi_noise(0.0, 0.02);
        for (const auto& o : obs) {
            const Eigen::Vector3d q = om_true.miller_to_lab(o.h, o.k, o.l);
            const double phi = nearest_ewald_phi_deg(q, lam, o.phi_deg);
            auto c = o;
            c.phi_obs_deg = phi + phi_noise(rng);
            obs_phi.push_back(c);
        }
        Eigen::Matrix3d rx;
        const double a = 0.2 * 3.14159265358979323846 / 180.0;
        rx << 1, 0, 0, 0, std::cos(a), -std::sin(a), 0, std::sin(a), std::cos(a);
        const auto om_x = OrientationMatrix::from_cell_and_rotation(cell0, rx * U);
        const auto rp = PositionalRefiner{}.refine(
            om_x, p4, lam, det0, bx_true, by_true, obs_phi
        );
        // Recovered orientation predicts the true φ for a sample reflection.
        const auto& o = obs_phi[100];
        const double phi_true = nearest_ewald_phi_deg(om_true.miller_to_lab(o.h, o.k, o.l), lam, o.phi_deg);
        const double phi_fit = nearest_ewald_phi_deg(rp.orientation.miller_to_lab(o.h, o.k, o.l), lam, o.phi_deg);
        std::cout << "phi-aware rms " << rp.rms_initial_px << " -> " << rp.rms_px
                  << " dphi=" << phi_fit - phi_true << '\n';
        CHECK(rp.ok && rp.improved);
        CHECK(std::fabs(phi_fit - phi_true) < 0.02);
        CHECK(std::fabs(rp.cell.c - 150.27) < 0.1);
    }

    // Too few observations → not ok, geometry unchanged.
    const std::vector<PositionalObservation> few(obs.begin(), obs.begin() + 10);
    const auto none = PositionalRefiner{}.refine(om0, p4, lam, det0, bx_true, by_true, few);
    CHECK(!none.ok && none.distance_mm == 200.0);

    return mosflm_test::finish("positional_refinement");
}
