#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "test_util.hpp"

#include <filesystem>
#include <vector>

using mosflm::diffraction::PredictedReflection;
using mosflm::diffraction::ReflectionPredictor;

namespace {

void compare(const std::filesystem::path& oracle,
             const std::vector<PredictedReflection>& got,
             const char* label) {
    CHECK_MSG(std::filesystem::exists(oracle),
              std::string(label) + ": missing " + oracle.string());
    const auto expected = ReflectionPredictor::read_dump(oracle);
    CHECK_EQ(got.size(), expected.size());
    double max_dx = 0.0;
    double max_dy = 0.0;
    const std::size_t n = std::min(got.size(), expected.size());
    for (std::size_t i = 0; i < n; ++i) {
        CHECK_EQ(got[i].h, expected[i].h);
        CHECK_EQ(got[i].k, expected[i].k);
        CHECK_EQ(got[i].l, expected[i].l);
        CHECK_EQ(got[i].full, expected[i].full);
        max_dx = std::max(max_dx, std::fabs(got[i].x_mm - expected[i].x_mm));
        max_dy = std::max(max_dy, std::fabs(got[i].y_mm - expected[i].y_mm));
    }
    CHECK_MSG(max_dx < 1e-6, std::string(label) + " max dx");
    CHECK_MSG(max_dy < 1e-6, std::string(label) + " max dy");
}

}  // namespace

int main() {
    const std::filesystem::path dir = "testdata/oracle/known_indexed";
    CHECK(std::filesystem::exists(dir));

    const auto still_result = ReflectionPredictor{}.predict(
        ReflectionPredictor::analytical_still_fixture()
    );
    const auto& still = still_result.reflections;
    compare(dir / "predict_still.dump", still, "known_indexed still");

    const auto osc_result = ReflectionPredictor{}.predict(
        ReflectionPredictor::analytical_oscillation_fixture()
    );
    const auto& osc = osc_result.reflections;
    compare(dir / "predict_osc.dump", osc, "known_indexed osc");

    return mosflm_test::finish("known_indexed_oracle");
}
