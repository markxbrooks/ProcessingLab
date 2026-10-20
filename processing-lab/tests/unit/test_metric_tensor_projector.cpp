/**
 * MetricTensorProjector: symmetry-derived unit-cell metric constraints.
 */
#include "mosflm/geometry/MetricTensorProjector.hpp"
#include "test_util.hpp"

#include <gemmi/symmetry.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

using mosflm::MetricTensorProjector;
using mosflm::UnitCell;

const gemmi::SpaceGroup* sg(int number) {
    const auto* s = gemmi::find_spacegroup_by_number(number);
    CHECK(s != nullptr);
    return s;
}

}  // namespace

int main() {
    using mosflm_test::close;

    // Free metric parameters follow from the rotations alone.
    struct Expect {
        int number;
        int free;
    };
    for (const Expect e : {Expect{1, 6}, Expect{3, 4}, Expect{5, 4},
                           Expect{16, 3}, Expect{75, 2}, Expect{89, 2},
                           Expect{79, 2}, Expect{143, 2}, Expect{146, 2},
                           Expect{168, 2}, Expect{195, 1}, Expect{196, 1}}) {
        const MetricTensorProjector p(sg(e.number));
        CHECK_MSG(
            p.free_parameters() == e.free,
            "SG " + std::to_string(e.number) + " free=" +
                std::to_string(p.free_parameters())
        );
    }

    // Critique example: P4 candidate, observed cell not exactly tetragonal.
    const UnitCell obs{57.810, 58.043, 150.016, 90.06, 90.07, 90.19};
    const MetricTensorProjector p4(sg(75));
    const auto r = p4.project(obs);
    CHECK(close(r.constrained.a, r.constrained.b, 1e-9));
    CHECK(close(r.constrained.alpha_deg, 90.0, 1e-9));
    CHECK(close(r.constrained.beta_deg, 90.0, 1e-9));
    CHECK(close(r.constrained.gamma_deg, 90.0, 1e-9));
    CHECK(close(r.constrained.c, 150.016, 1e-9));
    // Frobenius projection keeps a²+b²: RMS of a and b, ≈ arithmetic mean.
    const double rms = std::sqrt(0.5 * (obs.a * obs.a + obs.b * obs.b));
    CHECK(close(r.constrained.a, rms, 1e-9));
    CHECK(close(r.max_angle_dev_deg, 0.19, 1e-9));
    CHECK(r.dof == 4);
    CHECK(!r.chi2.has_value());
    std::cout << "P4 constrained a=" << r.constrained.a
              << " rel=" << r.rel_metric_dev << '\n';

    // Esd-weighted significance: 0.233 Å split is significant at 0.02 Å,
    // consistent at 0.3 Å.
    const mosflm::CellSigmas tight{0.02, 0.02, 0.05, 0.02, 0.02, 0.02};
    const auto rt = p4.project(obs, tight);
    CHECK(rt.chi2.has_value() && rt.p_value.has_value());
    CHECK_MSG(*rt.p_value < 1e-6, "tight esds should reject P4 metric");
    const mosflm::CellSigmas loose{0.3, 0.3, 0.5, 0.3, 0.3, 0.3};
    const auto rl = p4.project(obs, loose);
    CHECK_MSG(*rl.p_value > 0.05, "loose esds should accept P4 metric");
    // Unequal esds weight the better-determined axis.
    const mosflm::CellSigmas skew{0.01, 0.1, 0.05, 0.05, 0.05, 0.05};
    const auto rs = p4.project(obs, skew);
    CHECK(std::fabs(rs.constrained.a - obs.a) < std::fabs(rs.constrained.a - obs.b));

    // Hexagonal and rhombohedral (hexagonal axes) give γ = 120°.
    const UnitCell hex{75.1, 74.9, 120.0, 90.1, 89.9, 119.8};
    for (int n : {143, 146, 168}) {
        const auto h = MetricTensorProjector(sg(n)).project(hex);
        CHECK(close(h.constrained.a, h.constrained.b, 1e-9));
        CHECK(close(h.constrained.gamma_deg, 120.0, 1e-9));
        CHECK(close(h.constrained.alpha_deg, 90.0, 1e-9));
    }

    // Monoclinic unique-b: only α and γ fixed.
    const UnitCell mono{50.0, 60.0, 70.0, 90.2, 105.3, 89.7};
    const auto m = MetricTensorProjector(sg(5)).project(mono);
    CHECK(close(m.constrained.alpha_deg, 90.0, 1e-9));
    CHECK(close(m.constrained.gamma_deg, 90.0, 1e-9));
    CHECK(std::fabs(m.constrained.beta_deg - 105.3) < 0.05);

    // Cubic: a = b = c.
    const auto cu = MetricTensorProjector(sg(196)).project(
        UnitCell{100.1, 99.9, 100.3, 90.0, 90.1, 89.9}
    );
    CHECK(close(cu.constrained.a, cu.constrained.c, 1e-9));
    CHECK(close(cu.constrained.b, cu.constrained.c, 1e-9));

    // Triclinic is the identity.
    const auto tri = MetricTensorProjector(sg(1)).project(mono);
    CHECK(tri.max_edge_dev_A < 1e-9 && tri.max_angle_dev_deg < 1e-9);

    // Metric round trip and closed-form χ² tails.
    const UnitCell back = mosflm::cell_from_metric_tensor(mosflm::metric_tensor(mono));
    CHECK(close(back.beta_deg, mono.beta_deg, 1e-9));
    CHECK(close(mosflm::chi2_upper_tail(3.841458820694124, 1), 0.05, 1e-9));
    CHECK(close(mosflm::chi2_upper_tail(9.487729036781154, 4), 0.05, 1e-9));
    CHECK(close(mosflm::chi2_upper_tail(7.814727903251178, 3), 0.05, 1e-9));

    bool threw = false;
    try {
        (void)p4.project(obs, mosflm::CellSigmas{});
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    CHECK_MSG(threw, "zero esds must throw");

    return mosflm_test::finish("metric_tensor_projector");
}
