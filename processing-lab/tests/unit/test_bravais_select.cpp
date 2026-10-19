/**
 * SUGGLATT / LASTGOOD Bravais selection (penalty jump + highest crystal class).
 */

#include "mosflm/geometry/UnitCell.hpp"
#include "mosflm/indexing/BravaisClassifier.hpp"
#include "mosflm/indexing/BravaisSelect.hpp"
#include "test_util.hpp"

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

using mosflm::UnitCell;
using mosflm::indexing::BravaisCandidate;
using mosflm::indexing::BravaisClassifier;
using mosflm::indexing::crystal_class_rank;
using mosflm::indexing::find_lastgood_rank;
using mosflm::indexing::select_suggested_bravais;
using mosflm::indexing::symmetry_number_for_lattice;
using mosflm::indexing::write_bravais_penalty_table;

namespace {

BravaisCandidate make_cand(
    int character,
    const char* system_id,
    float criterion
) {
    BravaisCandidate c;
    c.character = character;
    c.system_id = system_id;
    c.criterion = criterion;
    c.orientation = Eigen::Matrix3d::Identity();
    return c;
}

}  // namespace

int main() {
    CHECK_EQ(crystal_class_rank("cI"), 6);
    CHECK_EQ(crystal_class_rank("tP"), 4);
    CHECK_EQ(crystal_class_rank("aP"), 1);
    CHECK_EQ(symmetry_number_for_lattice("tP"), 75);
    CHECK_EQ(symmetry_number_for_lattice("cI"), 197);
    CHECK_EQ(symmetry_number_for_lattice("aP"), 1);

    // Thaumatin-like reduced cell: a ≈ b ≠ c, angles 90°.
    const auto candidates =
        BravaisClassifier{}.classify_orthogonal_cell(57.98, 58.04, 148.16);
    CHECK_EQ(candidates.size(), 44u);

    // Old policy (min criterion) is almost always aP / character 44.
    int min_i = 0;
    float min_c = candidates[0].criterion;
    for (std::size_t i = 1; i < candidates.size(); ++i) {
        if (candidates[i].criterion < min_c) {
            min_c = candidates[i].criterion;
            min_i = static_cast<int>(i);
        }
    }
    CHECK_MSG(
        candidates[static_cast<std::size_t>(min_i)].system_id == "aP",
        "min-criterion should be aP"
    );

    // SUGGLATT / LASTGOOD: expect tetragonal (not aP).
    const std::size_t pick =
        select_suggested_bravais(candidates, 50.0f, 0, {});
    CHECK(pick < candidates.size());
    const auto& chosen = candidates[pick];
    CHECK_MSG(
        chosen.system_id == "tP" || chosen.system_id == "tI",
        "expected tetragonal SUGGLATT pick, got " + chosen.system_id +
            " crit=" + std::to_string(chosen.criterion)
    );
    CHECK(chosen.character != 44);
    const int sg = symmetry_number_for_lattice(chosen.system_id);
    CHECK(sg == 75 || sg == 79);  // tP → P4, tI → I4

    // prefer_system still wins.
    const std::size_t pref =
        select_suggested_bravais(candidates, 50.0f, 0, "aP");
    CHECK(pref < candidates.size());
    CHECK_MSG(candidates[pref].system_id == "aP", "prefer_system aP");

    // prefer_character 44.
    const std::size_t ch44 =
        select_suggested_bravais(candidates, 50.0f, 44, {});
    CHECK_EQ(static_cast<int>(ch44), 43);
    CHECK_EQ(candidates[ch44].character, 44);

    // Cubic cell → cubic family under SUGGLATT.
    const auto cubic =
        BravaisClassifier{}.classify_orthogonal_cell(78.19, 78.19, 78.19);
    const std::size_t cub =
        select_suggested_bravais(cubic, 50.0f, 0, {});
    CHECK(cub < cubic.size());
    CHECK_MSG(
        cubic[cub].system_id[0] == 'c',
        "expected cubic SUGGLATT pick, got " + cubic[cub].system_id
    );

    // Critique-like synthetic: low cluster 0..58 then jump to ~3279.
    // Soft cut 50 alone would drop tP@58; LASTGOOD must include it.
    std::vector<BravaisCandidate> critique = {
        make_cand(44, "aP", 0.0f),
        make_cand(43, "aP", 1.0f),
        make_cand(33, "mP", 7.0f),
        make_cand(34, "mP", 13.0f),
        make_cand(35, "mP", 45.0f),
        make_cand(32, "mC", 51.0f),
        make_cand(16, "oP", 51.0f),
        make_cand(17, "oC", 52.0f),
        make_cand(31, "mC", 53.0f),
        make_cand(11, "tP", 58.0f),
        make_cand(20, "mC", 3279.0f),
        make_cand(18, "oC", 3279.0f),
        make_cand(21, "mC", 3281.0f),
        make_cand(19, "oC", 3286.0f),
        make_cand(22, "mC", 3397.0f),
    };
    const int lastgood = find_lastgood_rank(critique);
    CHECK_MSG(
        lastgood == 10,
        "expected LASTGOOD=10 (before jump), got " + std::to_string(lastgood)
    );
    const std::size_t jump_pick =
        select_suggested_bravais(critique, 50.0f, 0, {});
    CHECK(jump_pick < critique.size());
    CHECK_MSG(
        critique[jump_pick].system_id == "tP",
        "critique jump should pick tP, got " + critique[jump_pick].system_id +
            " crit=" + std::to_string(critique[jump_pick].criterion)
    );
    CHECK(critique[jump_pick].criterion > 50.0f);

    // Mosflm-style penalty table header + LASTGOOD footer.
    UnitCell cell;
    cell.a = 57.98;
    cell.b = 58.04;
    cell.c = 148.16;
    cell.alpha_deg = 90.0;
    cell.beta_deg = 90.0;
    cell.gamma_deg = 90.0;
    std::ostringstream table;
    write_bravais_penalty_table(
        table,
        critique,
        critique[jump_pick].character,
        critique[jump_pick].system_id,
        cell,
        50.0f
    );
    const std::string text = table.str();
    CHECK_MSG(
        text.find("PENALTY LATT") != std::string::npos,
        "table missing PENALTY header"
    );
    CHECK_MSG(
        text.find("DPS Niggli goof") != std::string::npos,
        "table missing DPS goof note"
    );
    CHECK_MSG(
        text.find("LASTGOOD") != std::string::npos,
        "table missing LASTGOOD"
    );
    CHECK_MSG(
        text.find("Suggested Solution:") != std::string::npos,
        "table missing Suggested Solution"
    );
    CHECK_MSG(text.find("tP") != std::string::npos, "table missing tP");

    // Real-cell table still prints Suggested Solution.
    std::ostringstream table2;
    write_bravais_penalty_table(
        table2,
        candidates,
        chosen.character,
        chosen.system_id,
        cell,
        50.0f
    );
    CHECK_MSG(
        table2.str().find("Suggested Solution:") != std::string::npos,
        "thaumatin table missing Suggested Solution"
    );

    return mosflm_test::finish("bravais_select");
}
