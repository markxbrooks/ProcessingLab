/**
 * sum_partial_slices: rocking-curve slices → one reflection.
 */
#include "mosflm/diffraction/PartialSummation.hpp"
#include "test_util.hpp"

#include <cmath>

int main() {
    using mosflm::crystal::MosflmReflection;
    using mosflm::diffraction::sum_partial_slices;
    using mosflm_test::close;

    auto slice = [](int batch, float frac, float I, float phi, int flag) {
        MosflmReflection r;
        r.h = r.h0 = 1;
        r.k = r.k0 = 2;
        r.l = r.l0 = 3;
        r.batch = batch;
        r.rot = 10.0f;
        r.fractioncalc = frac;
        r.intensity = r.ipr = I;
        r.sigma = r.sigipr = 2.0f;
        r.phi_obs = phi;
        r.flag = flag;
        r.xdet_obs = 100.0f + static_cast<float>(batch);
        r.ydet_obs = 200.0f;
        return r;
    };
    std::vector<MosflmReflection> rows{
        slice(66, 0.25f, 25.0f, 9.925f, 0), slice(67, 0.5f, 50.0f, 10.075f, 0),
        slice(68, 0.24f, 24.0f, 10.225f, 4),  // weak-tail neg-profile bit dropped
        slice(68, 0.05f, -1.0f, 10.225f, 16),  // tiny tail: overlap bit ignored
    };
    // Same HKL 180° away is a different reflection.
    auto other = slice(1266, 1.0f, 7.0f, 190.0f, 0);
    other.rot = 190.0f;
    rows.push_back(other);

    const auto out = sum_partial_slices(rows);
    CHECK(out.size() == 2u);
    const auto& r = out[0];
    CHECK(r.batch == 67);
    CHECK(close(r.intensity, 98.0, 1e-4));
    CHECK(close(r.sigma, 4.0, 1e-4));
    CHECK(close(r.fractioncalc, 1.0, 1e-6));  // 1.04 capped
    CHECK(r.mpart == 0);
    CHECK(r.flag == 0);
    const double phi = (25 * 9.925 + 50 * 10.075 + 24 * 10.225) / 99.0;
    CHECK(close(r.phi_obs, phi, 1e-4));
    CHECK(close(r.xdet_obs, (25 * 166 + 50 * 167 + 24 * 168) / 99.0, 1e-3));
    CHECK(out[1].batch == 1266 && close(out[1].intensity, 7.0, 1e-6));

    // Incomplete edge reflection stays partial.
    const auto edge = sum_partial_slices({slice(1, 0.3f, 30.0f, 0.075f, 0)});
    CHECK(edge.size() == 1u && edge[0].mpart == 1 && close(edge[0].fractioncalc, 0.3, 1e-6));
    // Only a tail slice survived (core frames masked) → reflection dropped.
    const auto tail_only = sum_partial_slices({slice(70, 0.0005f, 0.0f, 10.5f, 0)});
    CHECK(tail_only.empty());
    return mosflm_test::finish("partial_summation");
}
