/**
 * IDXREF-v2: local indexing + geometry LS + INDEX_ORIGIN + distance refine.
 */

#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/geometry/LabFrame.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/indexing/BravaisSelect.hpp"
#include "mosflm/indexing/IdxrefRefine.hpp"
#include "mosflm/indexing/LocalIndexAssign.hpp"
#include "mosflm/indexing/SpotList.hpp"
#include "mosflm/refinement/PostRefinement.hpp"
#include "test_util.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

using mosflm::OrientationMatrix;
using mosflm::UnitCell;
using mosflm::indexing::IdxrefGeometry;
using mosflm::indexing::IdxrefObservation;
using mosflm::indexing::IdxrefSettings;
using mosflm::indexing::LocalIndexSettings;
using mosflm::indexing::Spot;
using mosflm::indexing::SpotList;
using mosflm::indexing::assign_local_indices;
using mosflm::indexing::build_idxref_observations;
using mosflm::indexing::lattice_metric_acceptable;
using mosflm::indexing::refine_idxref;
using mosflm::indexing::select_index_origin;
using mosflm::cpp_ub_to_mosflm_amat;

namespace {

bool near(double a, double b, double tol) {
    return std::fabs(a - b) <= tol;
}

}  // namespace

int main() {
    // --- Local indexing on a clean reciprocal lattice ---
    {
        UnitCell cell;
        cell.a = 60.0;
        cell.b = 80.0;
        cell.c = 100.0;
        const OrientationMatrix orient =
            OrientationMatrix::from_cell_and_rotation(cell);
        const Eigen::Matrix3d amat =
            cpp_ub_to_mosflm_amat(orient.ub(), 1.0);

        std::mt19937 rng(7);
        std::uniform_int_distribution<int> hkl(-4, 4);
        std::vector<float> x;
        std::vector<float> y;
        std::vector<float> z;
        for (int n = 0; n < 200; ++n) {
            Eigen::Vector3d h(
                static_cast<double>(hkl(rng)),
                static_cast<double>(hkl(rng)),
                static_cast<double>(hkl(rng))
            );
            if (h.squaredNorm() < 1e-12) {
                continue;
            }
            const Eigen::Vector3d rm = amat * h;
            x.push_back(static_cast<float>(rm.x()));
            y.push_back(static_cast<float>(rm.y()));
            z.push_back(static_cast<float>(rm.z()));
        }
        CHECK(x.size() >= 50u);
        LocalIndexSettings ls;
        const auto loc = assign_local_indices(amat, x, y, z, ls);
        CHECK(loc.n_indexed >= 40u);
        std::cout << "local_index: indexed=" << loc.n_indexed << "/"
                  << x.size() << " largest_sub=" << loc.n_largest_subtree
                  << '\n';
    }

    // --- build_idxref_observations maps image mm → MOSFLM XDET ---
    {
        IdxrefGeometry geo;
        geo.distance_mm = 200.0;
        geo.wavelength_A = 1.0;
        geo.beam_x_mm = 160.0;
        geo.beam_y_mm = 160.0;
        geo.pixel_mm = 0.075;
        geo.yscale = 1.0;

        SpotList spots;
        Spot s;
        s.x_mm = 175.5;  // image-frame absolute mm
        s.y_mm = 148.2;
        s.phi_end_deg = 12.0;
        s.intensity = 100.0;
        s.sigma = 1.0;
        spots.spots.push_back(s);

        mosflm::indexing::LocalIndexResult local;
        local.spots.resize(1);
        local.spots[0].indexed = true;
        local.spots[0].h = 1;
        local.spots[0].k = 0;
        local.spots[0].l = 0;
        local.n_indexed = 1;

        const auto obs = build_idxref_observations(
            spots, {0}, local, geo, nullptr, nullptr, nullptr
        );
        CHECK(obs.size() == 1u);
        // Identity through cpp_mm_to_mosflm_xdet with GENDATA undo:
        // XDET = x_mm / pixel, YDET = y_mm / pixel.
        CHECK(near(obs[0].xdet_obs_px, s.x_mm / geo.pixel_mm, 1e-6));
        CHECK(near(obs[0].ydet_obs_px, s.y_mm / geo.pixel_mm, 1e-6));
        std::cout << "idxref_xdet: OK\n";
    }

    // --- Geometry refine pulls a consistent wrong cell toward truth ---
    {
        UnitCell truth;
        truth.a = 58.0;
        truth.b = 58.0;
        truth.c = 150.0;
        const OrientationMatrix orient_true =
            OrientationMatrix::from_cell_and_rotation(truth);

        UnitCell start = truth;
        start.a = 62.0;
        start.b = 54.0;
        start.c = 142.0;
        start.beta_deg = 93.0;
        const OrientationMatrix orient_start =
            OrientationMatrix::from_cell_and_rotation(start);
        const Eigen::Matrix3d amat_start =
            cpp_ub_to_mosflm_amat(orient_start.ub(), 1.0);

        IdxrefGeometry geo;
        geo.distance_mm = 200.0;
        geo.wavelength_A = 1.0;
        geo.beam_x_mm = 160.0;
        geo.beam_y_mm = 160.0;
        geo.pixel_mm = 0.075;

        std::vector<IdxrefObservation> obs;
        std::mt19937 rng(11);
        std::uniform_int_distribution<int> hkl(-5, 5);
        std::normal_distribution<double> noise_px(0.0, 0.2);
        mosflm::DetectorGeometry det;
        det.distance_mm = geo.distance_mm;
        det.wavelength_A = geo.wavelength_A;
        det.pixel_size_x_mm = geo.pixel_mm;
        det.pixel_size_y_mm = geo.pixel_mm;
        for (int n = 0; n < 160; ++n) {
            const int h = hkl(rng);
            const int k = hkl(rng);
            const int l = hkl(rng);
            if (h == 0 && k == 0 && l == 0) {
                continue;
            }
            const Eigen::Vector3d q0 = orient_true.miller_to_lab(h, k, l);
            const auto sols =
                mosflm::refinement::ewald_phi_solutions_deg(q0, 1.0);
            if (sols.empty()) {
                continue;
            }
            const double phi = sols[0];
            double xdet = 0.0;
            double ydet = 0.0;
            if (!mosflm::refinement::predict_xdet_ydet_px(
                    orient_true, h, k, l, phi, geo.wavelength_A, det,
                    geo.beam_x_mm, geo.beam_y_mm, xdet, ydet
                )) {
                continue;
            }
            IdxrefObservation o;
            o.h = h;
            o.k = k;
            o.l = l;
            o.phi_obs_deg = phi;
            o.xdet_obs_px = xdet + noise_px(rng);
            o.ydet_obs_px = ydet + noise_px(rng);
            o.weight = 1.0;
            obs.push_back(o);
        }
        CHECK(obs.size() >= 30u);

        IdxrefSettings settings;
        settings.max_iterations = 16;
        settings.fd_cell_len = 0.15;
        settings.fd_cell_ang = 0.15;
        // Narrow φ span → F held fixed.
        settings.refine_distance_min_phi_span_deg = 1e9;
        const auto refined =
            refine_idxref(start, amat_start, obs, geo, settings);
        CHECK(refined.ok);
        std::cout << "idxref_refine: accepted=" << refined.accepted
                  << " cell=" << refined.cell.a << " " << refined.cell.b
                  << " " << refined.cell.c
                  << " rms_xy=" << refined.rms_pos_initial << "→"
                  << refined.rms_pos << '\n';
        CHECK(refined.accepted);
        CHECK(refined.rms_pos < refined.rms_pos_initial);
        CHECK(std::fabs(refined.cell.a - 58.0) < std::fabs(start.a - 58.0));
        CHECK(std::fabs(refined.cell.c - 150.0) < std::fabs(start.c - 150.0));
    }

    // --- INDEX_ORIGIN recovers a constant HKL offset ---
    {
        UnitCell cell;
        cell.a = 58.0;
        cell.b = 58.0;
        cell.c = 150.0;
        const OrientationMatrix orient =
            OrientationMatrix::from_cell_and_rotation(cell);
        const Eigen::Matrix3d amat =
            cpp_ub_to_mosflm_amat(orient.ub(), 1.0);

        IdxrefGeometry geo;
        geo.distance_mm = 200.0;
        geo.wavelength_A = 1.0;
        geo.beam_x_mm = 160.0;
        geo.beam_y_mm = 160.0;
        geo.pixel_mm = 0.075;

        mosflm::DetectorGeometry det;
        det.distance_mm = geo.distance_mm;
        det.wavelength_A = geo.wavelength_A;
        det.pixel_size_x_mm = geo.pixel_mm;
        det.pixel_size_y_mm = geo.pixel_mm;

        std::vector<IdxrefObservation> obs;
        std::mt19937 rng(13);
        std::uniform_int_distribution<int> hkl(-4, 4);
        for (int n = 0; n < 80; ++n) {
            const int h = hkl(rng);
            const int k = hkl(rng);
            const int l = hkl(rng);
            if (h == 0 && k == 0 && l == 0) {
                continue;
            }
            const Eigen::Vector3d q0 = orient.miller_to_lab(h, k, l);
            const auto sols =
                mosflm::refinement::ewald_phi_solutions_deg(q0, 1.0);
            if (sols.empty()) {
                continue;
            }
            const double phi = sols[0];
            double xdet = 0.0;
            double ydet = 0.0;
            if (!mosflm::refinement::predict_xdet_ydet_px(
                    orient, h, k, l, phi, geo.wavelength_A, det,
                    geo.beam_x_mm, geo.beam_y_mm, xdet, ydet
                )) {
                continue;
            }
            IdxrefObservation o;
            // Misindex by constant offset; RLP stays at the true lattice point.
            const Eigen::Vector3d rm = amat * Eigen::Vector3d(h, k, l);
            o.h = h + 1;
            o.k = k;
            o.l = l - 1;
            o.phi_obs_deg = phi;
            o.xdet_obs_px = xdet;
            o.ydet_obs_px = ydet;
            o.rx = rm.x();
            o.ry = rm.y();
            o.rz = rm.z();
            o.weight = 1.0;
            obs.push_back(o);
        }
        CHECK(obs.size() >= 20u);

        IdxrefSettings settings;
        settings.index_origin_max = 1;
        const Eigen::Vector3i origin =
            select_index_origin(cell, amat, obs, geo, settings);
        std::cout << "idxref_origin: " << origin.transpose() << '\n';
        CHECK(origin.x() == -1);
        CHECK(origin.y() == 0);
        CHECK(origin.z() == 1);
    }

    // --- Distance refine moves F toward truth on a wide φ wedge ---
    {
        UnitCell cell;
        cell.a = 58.0;
        cell.b = 58.0;
        cell.c = 150.0;
        const OrientationMatrix orient =
            OrientationMatrix::from_cell_and_rotation(cell);
        const Eigen::Matrix3d amat =
            cpp_ub_to_mosflm_amat(orient.ub(), 1.0);

        constexpr double kFTrue = 200.0;
        IdxrefGeometry geo;
        geo.distance_mm = 210.0;  // start 5% high
        geo.wavelength_A = 1.0;
        geo.beam_x_mm = 160.0;
        geo.beam_y_mm = 160.0;
        geo.pixel_mm = 0.075;

        mosflm::DetectorGeometry det_true;
        det_true.distance_mm = kFTrue;
        det_true.wavelength_A = geo.wavelength_A;
        det_true.pixel_size_x_mm = geo.pixel_mm;
        det_true.pixel_size_y_mm = geo.pixel_mm;

        std::vector<IdxrefObservation> obs;
        std::mt19937 rng(17);
        std::uniform_int_distribution<int> hkl(-5, 5);
        // Spread φ over >20° so distance refine is enabled.
        std::uniform_real_distribution<double> phi_pick(0.0, 40.0);
        for (int n = 0; n < 200; ++n) {
            const int h = hkl(rng);
            const int k = hkl(rng);
            const int l = hkl(rng);
            if (h == 0 && k == 0 && l == 0) {
                continue;
            }
            const Eigen::Vector3d q0 = orient.miller_to_lab(h, k, l);
            const auto sols =
                mosflm::refinement::ewald_phi_solutions_deg(q0, 1.0);
            if (sols.empty()) {
                continue;
            }
            // Prefer a solution near a random φ in the wedge.
            const double target = phi_pick(rng);
            double phi = sols[0];
            double best = std::fabs(sols[0] - target);
            for (double s : sols) {
                const double d = std::fabs(s - target);
                if (d < best) {
                    best = d;
                    phi = s;
                }
            }
            double xdet = 0.0;
            double ydet = 0.0;
            if (!mosflm::refinement::predict_xdet_ydet_px(
                    orient, h, k, l, phi, geo.wavelength_A, det_true,
                    geo.beam_x_mm, geo.beam_y_mm, xdet, ydet
                )) {
                continue;
            }
            IdxrefObservation o;
            o.h = h;
            o.k = k;
            o.l = l;
            o.phi_obs_deg = phi;
            o.xdet_obs_px = xdet;
            o.ydet_obs_px = ydet;
            o.weight = 1.0;
            obs.push_back(o);
        }
        CHECK(obs.size() >= 40u);

        IdxrefSettings settings;
        settings.max_iterations = 20;
        settings.refine_distance_min_phi_span_deg = 20.0;
        settings.fd_distance_mm = 0.25;
        settings.distance_frac_bound = 0.08;
        // Hold cell fixed enough that F can move.
        settings.fd_cell_len = 0.02;
        settings.fd_cell_ang = 0.02;
        const auto refined =
            refine_idxref(cell, amat, obs, geo, settings);
        CHECK(refined.ok);
        std::cout << "idxref_distance: F=" << refined.distance_mm_initial
                  << "→" << refined.distance_mm
                  << " rms_xy=" << refined.rms_pos_initial << "→"
                  << refined.rms_pos << '\n';
        CHECK(refined.accepted);
        CHECK(std::fabs(refined.distance_mm - kFTrue) <
              std::fabs(geo.distance_mm - kFTrue));
    }

    // --- Bravais metric gate rejects skewed mC ---
    {
        UnitCell skewed;
        skewed.a = 55.0;
        skewed.b = 59.0;
        skewed.c = 163.0;
        skewed.alpha_deg = 44.0;
        skewed.beta_deg = 138.0;
        skewed.gamma_deg = 102.0;
        CHECK(!lattice_metric_acceptable("mC", skewed));
        CHECK(!lattice_metric_acceptable("tP", skewed));
        CHECK(lattice_metric_acceptable("aP", skewed));

        UnitCell tet;
        tet.a = 58.0;
        tet.b = 58.0;
        tet.c = 150.0;
        CHECK(lattice_metric_acceptable("tP", tet));
        std::cout << "idxref_bravais_gate: OK\n";
    }

    return mosflm_test::finish("idxref_refine");
}
