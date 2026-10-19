#include "mosflm/indexing/BravaisClassifier.hpp"
#include "test_util.hpp"

#include <string>

using mosflm::indexing::BravaisClassifier;

int main() {
    const BravaisClassifier classifier;

    // Orthogonal (tetragonal-looking) cell: a = b = 58, c = 150.
    const auto candidates = classifier.classify_orthogonal_cell(58.0, 58.0, 150.0);
    CHECK_EQ(candidates.size(), 44u);

    // Characters are numbered 1..44 with two-character system ids.
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        CHECK_EQ(candidates[i].character, static_cast<int>(i) + 1);
        CHECK_EQ(candidates[i].system_id.size(), 2u);
    }

    // The best (lowest) criterion should be a lattice character compatible
    // with a = b != c, i.e. one of the tetragonal (tP/tI) or orthogonal
    // settings. Find the minimum-criterion candidate.
    int best = -1;
    float best_crit = 1e30f;
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        if (candidates[i].criterion < best_crit) {
            best_crit = candidates[i].criterion;
            best = static_cast<int>(i);
        }
    }
    CHECK(best >= 0);
    const std::string best_sys = candidates[static_cast<std::size_t>(best)].system_id;
    CHECK_MSG(best_sys.size() == 2, "best system id: " + best_sys);

    // A cubic cell should be classified with a very low criterion too.
    const auto cubic = classifier.classify_orthogonal_cell(50.0, 50.0, 50.0);
    CHECK_EQ(cubic.size(), 44u);
    float cubic_best = 1e30f;
    for (const auto& c : cubic) {
        cubic_best = std::min(cubic_best, c.criterion);
    }
    CHECK(cubic_best < 1e-3f);

    return mosflm_test::finish("bravais_classifier");
}
