#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "test_util.hpp"

#include <cstdio>
#include <filesystem>
#include <vector>

using mosflm::diffraction::PredictedReflection;
using mosflm::diffraction::ReflectionPredictor;

int main() {
    const char* env = std::getenv("MOSFLM_PREDICT_OSC_ORACLE");
    if (env == nullptr) {
        std::cerr << "MOSFLM_PREDICT_OSC_ORACLE not set\n";
        return 1;
    }
    const std::filesystem::path oracle(env);
    const std::vector<PredictedReflection> expected =
        ReflectionPredictor::read_dump(oracle);

    const auto req = ReflectionPredictor::analytical_oscillation_fixture();
    const auto result = ReflectionPredictor{}.predict(req);
    const std::vector<PredictedReflection>& got = result.reflections;

    CHECK_EQ(got.size(), expected.size());
    CHECK(!got.empty());
    CHECK_EQ(expected.size(), 368u);

    std::size_t n_partial = 0;
    double max_dx = 0.0;
    const std::size_t n = std::min(got.size(), expected.size());
    for (std::size_t i = 0; i < n; ++i) {
        CHECK_EQ(got[i].h, expected[i].h);
        CHECK_EQ(got[i].k, expected[i].k);
        CHECK_EQ(got[i].l, expected[i].l);
        CHECK_EQ(got[i].full, expected[i].full);
        max_dx = std::max(max_dx, std::fabs(got[i].x_mm - expected[i].x_mm));
        if (!got[i].full) {
            ++n_partial;
        }
    }
    CHECK_CLOSE(max_dx, 0.0, 1e-6);
    CHECK_EQ(n_partial, 368u);

    // Oscillation window endpoints are within the requested scan.
    for (const auto& s : got) {
        CHECK(s.phi_begin_deg >= req.scan.phi_start_deg - 1e-9);
        CHECK(s.phi_end_deg <= req.scan.phi_end_deg + 1e-9);
    }

    return mosflm_test::finish("predict_osc");
}
