#include "mosflm/control/ControlSession.hpp"
#include "mosflm/diffraction/Strategy.hpp"
#include "mosflm/io/RunRecord.hpp"
#include "test_util.hpp"

#include <Eigen/Geometry>

#include <cmath>
#include <iostream>

namespace {

using mosflm::diffraction::PredictionRequest;
using mosflm::diffraction::StrategySettings;

/// A 50 Å cubic crystal in an arbitrary orientation, and a detector large
/// enough to record everything to 2.5 Å.
PredictionRequest geometry(const char* space_group) {
    const mosflm::UnitCell cell{50.0, 50.0, 50.0, 90.0, 90.0, 90.0};
    PredictionRequest req;
    req.crystal = mosflm::crystal::CrystalMetadata::from_mosflm_cell(cell, space_group);
    const Eigen::Matrix3d rot =
        Eigen::AngleAxisd(0.6, Eigen::Vector3d(0.2, 0.7, -0.4).normalized()).toRotationMatrix();
    req.orientation = mosflm::OrientationMatrix::from_cell_and_rotation(cell, rot);
    req.detector.wavelength_A = 1.0;
    req.detector.distance_mm = 100.0;
    req.detector.pixel_size_x_mm = 0.172;
    req.detector.pixel_size_y_mm = 0.172;
    req.detector.nfast = 3000;
    req.detector.nslow = 3000;
    req.settings.dmin_A = 2.5;
    req.settings.dmax_A = 100.0;
    return req;
}

}  // namespace

int main() {
    using mosflm::diffraction::compute_strategy;

    StrategySettings s;
    s.target_completeness = 0.90;
    s.osc_deg = 0.25;
    s.rotation_deg = 30.0;
    const auto p1 = compute_strategy(geometry("P 1"), s);
    CHECK(p1.unique_possible > 3000);
    CHECK_EQ(p1.unique_possible, p1.acentric_possible);  // P1: no centric reflections
    CHECK(p1.dmax_A < 60.0 && p1.dmin_A == 2.5);
    // The best completeness never falls as the rotation grows; anomalous
    // pairs need both mates, so lag behind.
    for (std::size_t i = 1; i < p1.curve.size(); ++i) {
        CHECK(p1.curve[i].completeness >= p1.curve[i - 1].completeness - 1e-12);
        CHECK(p1.curve[i].anomalous_completeness <= p1.curve[i].completeness + 1e-12);
    }
    CHECK(p1.target.has_value());
    if (p1.target) {
        // P1 needs most of a half turn, less the blind region's share.
        CHECK(p1.target->rotation_deg > 100.0 && p1.target->rotation_deg < 200.0);
        CHECK(p1.target->completeness >= 0.90);
        CHECK_EQ(p1.target->images, static_cast<int>(std::ceil(p1.target->rotation_deg / 0.25 - 1e-6)));
        CHECK(!p1.target->shells.empty() && p1.target->multiplicity >= 1.0);
        int possible = 0;
        for (const auto& sh : p1.target->shells) {
            possible += sh.possible;
        }
        CHECK_EQ(possible, p1.unique_possible);
    }
    CHECK(p1.rotation.has_value() && std::fabs(p1.rotation->rotation_deg - 30.0) < 1e-9);
    CHECK(p1.max_completeness > 0.95);

    // Cubic symmetry: the same completeness in a much smaller rotation.
    const auto cubic = compute_strategy(geometry("P 4 3 2"), s);
    CHECK(cubic.unique_possible < p1.unique_possible / 10);
    CHECK(cubic.target.has_value());
    if (cubic.target && p1.target) {
        CHECK(cubic.target->rotation_deg < p1.target->rotation_deg / 3.0);
    }

    // Anomalous: Friedel pairs need more rotation.
    StrategySettings anom = s;
    anom.anomalous = true;
    const auto p1a = compute_strategy(geometry("P 1"), anom);
    if (p1a.target && p1.target) {
        CHECK(p1a.target->rotation_deg > p1.target->rotation_deg);
        CHECK(p1a.target->anomalous_completeness >= 0.90);
    } else {
        CHECK(p1a.max_anomalous_completeness < 0.90);
    }

    // Within a collected 60° sweep: windows stay inside it.
    StrategySettings sweep = s;
    sweep.sweep_deg = std::pair(10.0, 70.0);
    sweep.target_completeness = 0.30;
    const auto in = compute_strategy(geometry("P 1"), sweep);
    CHECK(!in.curve.empty() && in.curve.back().rotation_deg <= 60.0 + 1e-9);
    CHECK(in.target.has_value());
    if (in.target) {
        CHECK(in.target->start_deg >= 10.0 - 1e-9 && in.target->end_deg <= 70.0 + 1e-9);
    }

    // CONTROL STRATEGY keywords.
    {
        mosflm::control::ControlSession session;
        session.parse_line("STRATEGY AUTO ROTATE 60 COMPLETENESS 95 ANOM FRAMES 900 1 RESOLUTION 2.1 STEP 2");
        const auto& st = session.state();
        CHECK(st.strategy_pending && st.strategy.anomalous);
        CHECK_CLOSE(st.strategy.rotation_deg, 60.0, 1e-12);
        CHECK_CLOSE(st.strategy.target_completeness, 0.95, 1e-12);
        CHECK(st.strategy_frames == std::pair(1, 900));
        CHECK_CLOSE(st.strategy_dmin_A, 2.1, 1e-12);
        CHECK_CLOSE(st.strategy.step_deg, 2.0, 1e-12);
        session.parse_line("STRATEGY FULL NOANOM");
        CHECK(!session.state().strategy_frames && !session.state().strategy.anomalous);
        bool threw = false;
        try {
            session.parse_line("STRATEGY SIDEWAYS");
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        CHECK(threw);
        CHECK(mosflm::io::planned_stages("XPARM x\nSTRATEGY COMPLETENESS 99\nGO\nEXIT\n") ==
              std::vector<mosflm::io::RunStep>{mosflm::io::RunStep::STRATEGY});
    }

    return mosflm_test::finish("strategy");
}
