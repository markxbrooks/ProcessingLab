#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "test_util.hpp"

#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

using mosflm::diffraction::PredictedReflection;
using mosflm::diffraction::ReflectionPredictor;

int main() {
    const char* env = std::getenv("MOSFLM_PREDICT_ORACLE");
    if (env == nullptr) {
        std::cerr << "MOSFLM_PREDICT_ORACLE not set\n";
        return 1;
    }
    const std::filesystem::path oracle(env);
    std::vector<PredictedReflection> expected =
        ReflectionPredictor::read_dump(oracle);

    const auto req = ReflectionPredictor::analytical_still_fixture();
    const auto result = ReflectionPredictor{}.predict(req);
    const std::vector<PredictedReflection>& got = result.reflections;

    CHECK_EQ(got.size(), expected.size());
    CHECK(!got.empty());

    std::size_t n = std::min(got.size(), expected.size());
    double max_dx = 0.0;
    double max_dy = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        CHECK_EQ(got[i].h, expected[i].h);
        CHECK_EQ(got[i].k, expected[i].k);
        CHECK_EQ(got[i].l, expected[i].l);
        CHECK_EQ(got[i].full, expected[i].full);
        max_dx = std::max(max_dx, std::fabs(got[i].x_mm - expected[i].x_mm));
        max_dy = std::max(max_dy, std::fabs(got[i].y_mm - expected[i].y_mm));
    }
    CHECK_MSG(max_dx < 1e-6, "max dx = " + std::to_string(max_dx));
    CHECK_MSG(max_dy < 1e-6, "max dy = " + std::to_string(max_dy));

    // write_dump / read_dump roundtrip.
    const std::filesystem::path tmp =
        std::filesystem::temp_directory_path() / "mosflm_predict_stub.dump";
    ReflectionPredictor::write_dump(tmp, got);
    const auto reread = ReflectionPredictor::read_dump(tmp);
    CHECK_EQ(reread.size(), got.size());
    if (!reread.empty()) {
        CHECK_EQ(reread[0].h, got[0].h);
        CHECK_CLOSE(reread[0].x_mm, got[0].x_mm, 1e-5);
    }
    std::remove(tmp.c_str());

    // All stills predictions are "full".
    for (const auto& s : got) {
        CHECK(s.full);
    }

    return mosflm_test::finish("predict_stub");
}
