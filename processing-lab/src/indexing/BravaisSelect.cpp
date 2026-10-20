#include "mosflm/indexing/BravaisSelect.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <vector>

extern "C" {
#include "util.h"
#include "cryst_matrix.h"
}

namespace mosflm::indexing {
namespace {

UnitCell cell_from_dps_om(const Eigen::Matrix3d& orientation) {
    struct matrix_3x3 om{};
    for (int i = 0; i < 3; ++i) {
        om.r[i].x = static_cast<float>(orientation(i, 0));
        om.r[i].y = static_cast<float>(orientation(i, 1));
        om.r[i].z = static_cast<float>(orientation(i, 2));
    }
    struct cell_struct cell{};
    dircell_from_om(om, &cell);
    UnitCell out;
    out.a = cell.a;
    out.b = cell.b;
    out.c = cell.c;
    constexpr double kRad = 180.0 / 3.14159265358979323846;
    out.alpha_deg = cell.alpha * kRad;
    out.beta_deg = cell.beta * kRad;
    out.gamma_deg = cell.gamma * kRad;
    return out;
}

/// Mosflm LTYPE / SPLIST / SPGMIN parallel tables (lattice index 0..14).
struct LatticeLabels {
    const char* system;
    const char* spgmin;
    const char* possible;
};

constexpr LatticeLabels kLatticeLabels[] = {
    {"aP", "P1", "P1"},
    {"mP", "P2", "P2,P21"},
    {"mC", "C2", "C2"},
    {"mI", "C2", "C2 (transformed from I2)"},
    {"oP", "P222", "P222,P2221,P21212,P212121"},
    {"oC", "C222", "C222,C2221"},
    {"oF", "F222", "F222"},
    {"oI", "I222", "I222,I212121"},
    {"tP",
     "P4",
     "P4,P41,P42,P43,P422,P4212,P4122,P41212,P4222,P42212,P4322,P43212"},
    {"tI", "I4", "I4,I41,I422,I4122"},
    {"hP",
     "P3",
     "P3,P31,P32,P312,P321,P3112,P3121,P3212,P3221"},
    {"hR", "H3", "H3,H32 (hexagonal settings of R3 and R32)"},
    {"cP", "P23", "P23,P213,P432,P4232,P4332,P4132"},
    {"cF", "F23", "F23,F432,F4132"},
    {"cI", "I23", "I23,I213,I432,I4132"},
};

constexpr const char* kHexagonalP6 =
    "P6,P61,P65,P62,P64,P63,P622,P6122,P6522,P6222,P6422,P6322";

const LatticeLabels* find_lattice(std::string_view system_id) {
    for (const auto& row : kLatticeLabels) {
        if (system_id == row.system) {
            return &row;
        }
    }
    return nullptr;
}

/// Indices sorted by ascending DPS goof (Mosflm IRESORDER); tie-break character.
std::vector<std::size_t> order_by_penalty(
    const std::vector<BravaisCandidate>& candidates
) {
    std::vector<std::size_t> order(candidates.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
    }
    std::stable_sort(
        order.begin(),
        order.end(),
        [&](std::size_t a, std::size_t b) {
            if (candidates[a].criterion != candidates[b].criterion) {
                return candidates[a].criterion < candidates[b].criterion;
            }
            return candidates[a].character < candidates[b].character;
        }
    );
    return order;
}

/// Mosflm MERIT scale for jump detection only (printed PENALTY stays raw).
int merit_from_criterion(float criterion, float merit_max) {
    if (!(merit_max > 0.0f)) {
        return 0;
    }
    return static_cast<int>(
        std::abs(std::lround(static_cast<double>(criterion) * 999.0 /
                             static_cast<double>(merit_max)))
    );
}

float merit_max_of(const std::vector<BravaisCandidate>& candidates) {
    float m = 0.0f;
    for (const auto& c : candidates) {
        if (c.criterion > m) {
            m = c.criterion;
        }
    }
    return m;
}

/// Highest crystal-class among sorted ranks ``1..lastgood`` (1-based inclusive).
std::size_t pick_highest_symmetry_in_window(
    const std::vector<BravaisCandidate>& candidates,
    const std::vector<std::size_t>& order,
    int lastgood
) {
    if (order.empty() || lastgood < 1) {
        return static_cast<std::size_t>(-1);
    }
    const int lim = std::min(lastgood, static_cast<int>(order.size()));
    std::size_t best = static_cast<std::size_t>(-1);
    int best_rank = -1;
    float best_crit = std::numeric_limits<float>::max();
    for (int rank = 1; rank <= lim; ++rank) {
        const std::size_t i = order[static_cast<std::size_t>(rank - 1)];
        const auto& c = candidates[i];
        const int cr = crystal_class_rank(c.system_id);
        if (cr > best_rank ||
            (cr == best_rank && c.criterion < best_crit) ||
            (cr == best_rank && c.criterion == best_crit &&
             (best == static_cast<std::size_t>(-1) ||
              c.character < candidates[best].character))) {
            best_rank = cr;
            best_crit = c.criterion;
            best = i;
        }
    }
    return best;
}

std::size_t pick_soft_cut(
    const std::vector<BravaisCandidate>& candidates,
    float max_penalty
) {
    std::size_t best = static_cast<std::size_t>(-1);
    int best_rank = -1;
    float best_crit = std::numeric_limits<float>::max();
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        const auto& c = candidates[i];
        if (c.criterion > max_penalty) {
            continue;
        }
        const int rank = crystal_class_rank(c.system_id);
        if (rank > best_rank ||
            (rank == best_rank && c.criterion < best_crit) ||
            (rank == best_rank && c.criterion == best_crit &&
             (best == static_cast<std::size_t>(-1) ||
              c.character < candidates[best].character))) {
            best_rank = rank;
            best_crit = c.criterion;
            best = i;
        }
    }
    return best;
}

std::size_t pick_min_criterion(
    const std::vector<BravaisCandidate>& candidates
) {
    std::size_t best = 0;
    float best_crit = candidates[0].criterion;
    for (std::size_t i = 1; i < candidates.size(); ++i) {
        if (candidates[i].criterion < best_crit) {
            best_crit = candidates[i].criterion;
            best = i;
        }
    }
    return best;
}

}  // namespace

int crystal_class_rank(std::string_view system_id) {
    if (system_id.empty()) {
        return 0;
    }
    const char c = static_cast<char>(
        std::tolower(static_cast<unsigned char>(system_id[0]))
    );
    switch (c) {
        case 'c':
            return 6;  // cubic
        case 'h':
        case 'r':
            return 5;  // hexagonal / rhombohedral
        case 't':
            return 4;  // tetragonal
        case 'o':
            return 3;  // orthorhombic
        case 'm':
            return 2;  // monoclinic
        case 'a':
            return 1;  // triclinic (aP)
        default:
            return 0;
    }
}

int symmetry_number_for_lattice(std::string_view system_id) {
    // Mosflm SPGMIN / DiffractionLab DEFAULT_SPACEGROUP_FOR_LATTICE.
    if (system_id == "aP") {
        return 1;
    }
    if (system_id == "mP") {
        return 3;
    }
    if (system_id == "mC" || system_id == "mI") {
        return 5;
    }
    if (system_id == "oP") {
        return 16;
    }
    if (system_id == "oC") {
        return 21;
    }
    if (system_id == "oF") {
        return 22;
    }
    if (system_id == "oI") {
        return 23;
    }
    if (system_id == "tP") {
        return 75;
    }
    if (system_id == "tI") {
        return 79;
    }
    if (system_id == "hP") {
        return 143;
    }
    if (system_id == "hR") {
        return 146;
    }
    if (system_id == "cP") {
        return 195;
    }
    if (system_id == "cF") {
        return 196;
    }
    if (system_id == "cI") {
        return 197;
    }
    return 1;
}

std::string_view spgmin_label_for_lattice(std::string_view system_id) {
    if (const LatticeLabels* row = find_lattice(system_id)) {
        return row->spgmin;
    }
    return "P1";
}

std::string_view possible_spacegroups_for_lattice(std::string_view system_id) {
    if (const LatticeLabels* row = find_lattice(system_id)) {
        return row->possible;
    }
    return "P1";
}

int find_lastgood_rank(const std::vector<BravaisCandidate>& candidates) {
    // Mosflm DPSCHOOSE: LASTGOOD starts at 1; scan ranks 2..N-1 for the
    // largest relative MERIT jump with MERIT(i) ≤ 100 (penalty-only DOF=1).
    if (candidates.size() < 2) {
        return 1;
    }
    const auto order = order_by_penalty(candidates);
    const int n = static_cast<int>(order.size());
    const float merit_max = merit_max_of(candidates);

    int lastgood = 1;
    float maxgood = 0.1f;
    for (int i = 2; i <= n - 1; ++i) {
        const std::size_t idx = order[static_cast<std::size_t>(i - 1)];
        const int merit_i = merit_from_criterion(
            candidates[idx].criterion, merit_max
        );
        float testgood = 0.0f;
        if (merit_i <= 100) {
            // Second solution is often another aP with MERIT 0 — ignore.
            if (i == 2 && merit_i == 0 &&
                candidates[idx].system_id == "aP") {
                testgood = maxgood;
            } else {
                const std::size_t idx_next =
                    order[static_cast<std::size_t>(i)];
                const int merit_next = merit_from_criterion(
                    candidates[idx_next].criterion, merit_max
                );
                const int next_cap = std::min(100, merit_next);
                testgood = (static_cast<float>(next_cap) -
                            static_cast<float>(merit_i)) /
                           (static_cast<float>(merit_i) + 1.0f);
            }
        }
        if (testgood > maxgood) {
            maxgood = testgood;
            lastgood = i;
        }
    }
    return lastgood;
}

namespace {

bool near_angle(double value_deg, double ideal_deg, double tol_deg) {
    return std::fabs(value_deg - ideal_deg) <= tol_deg;
}

bool near_axis_pair(double a, double b, double max_frac) {
    const double m = 0.5 * (std::fabs(a) + std::fabs(b));
    if (m < 1e-6) {
        return false;
    }
    return std::fabs(a - b) <= max_frac * m;
}

}  // namespace

bool lattice_metric_acceptable(
    std::string_view system_id,
    const UnitCell& cell,
    double max_axis_frac,
    double max_angle_deg
) {
    if (!cell.is_valid() || system_id.empty()) {
        return false;
    }
    if (system_id == "aP") {
        return true;
    }

    const double aa = cell.alpha_deg;
    const double bb = cell.beta_deg;
    const double gg = cell.gamma_deg;
    const bool orth =
        near_angle(aa, 90.0, max_angle_deg) &&
        near_angle(bb, 90.0, max_angle_deg) &&
        near_angle(gg, 90.0, max_angle_deg);

    if (system_id.size() < 2) {
        return false;
    }
    const char fam = static_cast<char>(
        std::tolower(static_cast<unsigned char>(system_id[0]))
    );
    const char cent = static_cast<char>(
        std::tolower(static_cast<unsigned char>(system_id[1]))
    );

    if (fam == 'm') {
        // Monoclinic: unique angle (β) free; α=γ=90°.
        return near_angle(aa, 90.0, max_angle_deg) &&
               near_angle(gg, 90.0, max_angle_deg);
    }
    if (fam == 'o') {
        return orth;
    }
    if (fam == 't') {
        return orth && near_axis_pair(cell.a, cell.b, max_axis_frac);
    }
    if (fam == 'h' || fam == 'r') {
        // Hexagonal setting: a≈b, α=β=90°, γ=120°.
        if (cent == 'r') {
            // Rhombohedral: a≈b≈c, angles equal (soft).
            return near_axis_pair(cell.a, cell.b, max_axis_frac) &&
                   near_axis_pair(cell.b, cell.c, max_axis_frac) &&
                   near_angle(aa, bb, max_angle_deg) &&
                   near_angle(bb, gg, max_angle_deg);
        }
        return near_axis_pair(cell.a, cell.b, max_axis_frac) &&
               near_angle(aa, 90.0, max_angle_deg) &&
               near_angle(bb, 90.0, max_angle_deg) &&
               near_angle(gg, 120.0, max_angle_deg);
    }
    if (fam == 'c') {
        return orth && near_axis_pair(cell.a, cell.b, max_axis_frac) &&
               near_axis_pair(cell.b, cell.c, max_axis_frac);
    }
    return false;
}

std::size_t select_bravais_after_idxref(
    const std::vector<BravaisCandidate>& candidates,
    const UnitCell& refined_cell,
    float max_penalty,
    int prefer_character,
    std::string_view prefer_system
) {
    if (candidates.empty()) {
        return static_cast<std::size_t>(-1);
    }

    // Explicit user preference (insulin SOLUTION / SYSTEM) wins when present.
    if (prefer_character >= 1 && prefer_character <= 44) {
        const std::size_t idx = static_cast<std::size_t>(prefer_character - 1);
        if (idx < candidates.size()) {
            return idx;
        }
    }
    if (!prefer_system.empty()) {
        float best_c = std::numeric_limits<float>::max();
        std::size_t best = static_cast<std::size_t>(-1);
        for (std::size_t i = 0; i < candidates.size(); ++i) {
            if (candidates[i].system_id == prefer_system &&
                candidates[i].criterion < best_c) {
                best_c = candidates[i].criterion;
                best = i;
            }
        }
        if (best != static_cast<std::size_t>(-1)) {
            return best;
        }
    }

    // Default: stay on refined triclinic aP (XDS continues aP until CORRECT).
    std::size_t ap_idx = static_cast<std::size_t>(-1);
    float ap_crit = std::numeric_limits<float>::max();
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        if (candidates[i].system_id == "aP" &&
            candidates[i].criterion < ap_crit) {
            ap_crit = candidates[i].criterion;
            ap_idx = i;
        }
    }

    // Accept higher symmetry only if metric matches and criterion is not junk.
    float crit_cut =
        std::max(max_penalty, ap_crit < 1e20f ? (ap_crit + 25.0f) : 25.0f);
    // A Mosflm LASTGOOD penalty jump admits everything below it, as in
    // select_suggested_bravais (real-data raw goofs sit far above 50).
    const int lastgood = find_lastgood_rank(candidates);
    if (lastgood > 1) {
        const auto order = order_by_penalty(candidates);
        crit_cut = std::max(
            crit_cut,
            candidates[order[static_cast<std::size_t>(lastgood - 1)]].criterion
        );
    }
    std::size_t best = ap_idx;
    int best_rank = (ap_idx != static_cast<std::size_t>(-1))
                        ? crystal_class_rank("aP")
                        : 0;
    float best_crit = (ap_idx != static_cast<std::size_t>(-1))
                          ? candidates[ap_idx].criterion
                          : std::numeric_limits<float>::max();

    for (std::size_t i = 0; i < candidates.size(); ++i) {
        const auto& c = candidates[i];
        if (c.system_id == "aP") {
            continue;
        }
        if (c.criterion > crit_cut) {
            continue;
        }
        // Gate on the candidate's own conventional cell: the aP cell's axis
        // order (and centring) need not match the candidate setting.
        const UnitCell own = cell_from_dps_om(c.orientation);
        if (!lattice_metric_acceptable(
                c.system_id, own.is_valid() ? own : refined_cell
            )) {
            continue;
        }
        const int rank = crystal_class_rank(c.system_id);
        // Prefer t/o/h/c over m: higher crystal_class_rank wins; tie → lower crit.
        if (rank > best_rank ||
            (rank == best_rank && c.criterion < best_crit)) {
            best_rank = rank;
            best_crit = c.criterion;
            best = i;
        }
    }

    if (best != static_cast<std::size_t>(-1)) {
        return best;
    }
    return pick_min_criterion(candidates);
}

std::size_t select_suggested_bravais(
    const std::vector<BravaisCandidate>& candidates,
    float max_penalty,
    int prefer_character,
    std::string_view prefer_system
) {
    if (candidates.empty()) {
        return static_cast<std::size_t>(-1);
    }

    if (prefer_character >= 1 && prefer_character <= 44) {
        const std::size_t idx = static_cast<std::size_t>(prefer_character - 1);
        if (idx < candidates.size()) {
            return idx;
        }
    }

    if (!prefer_system.empty()) {
        float best_c = std::numeric_limits<float>::max();
        std::size_t best = static_cast<std::size_t>(-1);
        for (std::size_t i = 0; i < candidates.size(); ++i) {
            if (candidates[i].system_id == prefer_system &&
                candidates[i].criterion < best_c) {
                best_c = candidates[i].criterion;
                best = i;
            }
        }
        if (best != static_cast<std::size_t>(-1)) {
            return best;
        }
    }

    // Mosflm LASTGOOD jump window (when a discontinuity was found).
    const int lastgood = find_lastgood_rank(candidates);
    if (lastgood > 1) {
        const auto order = order_by_penalty(candidates);
        const std::size_t jump_pick =
            pick_highest_symmetry_in_window(candidates, order, lastgood);
        if (jump_pick != static_cast<std::size_t>(-1)) {
            return jump_pick;
        }
    }

    // Soft-cut fallback: highest crystal class among penalty ≤ max_penalty.
    const std::size_t soft = pick_soft_cut(candidates, max_penalty);
    if (soft != static_cast<std::size_t>(-1)) {
        return soft;
    }

    // Global minimum criterion (typically aP / 44).
    return pick_min_criterion(candidates);
}

void write_bravais_penalty_table(
    std::ostream& os,
    const std::vector<BravaisCandidate>& candidates,
    int selected_character,
    std::string_view selected_system,
    const UnitCell& selected_cell,
    float max_penalty
) {
    if (candidates.empty()) {
        return;
    }

    const auto order = order_by_penalty(candidates);
    const int lastgood = find_lastgood_rank(candidates);

    int suggested_rank = 0;
    float suggested_penalty = 0.0f;
    for (std::size_t rank = 0; rank < order.size(); ++rank) {
        const auto& c = candidates[order[rank]];
        if (c.character == selected_character) {
            suggested_rank = static_cast<int>(rank) + 1;
            suggested_penalty = c.criterion;
            break;
        }
    }

    os << " List of possible Laue groups, sorted on penalty index.\n"
       << " The lower the PENALTY, the better\n"
       << " PENALTY = DPS Niggli goof (raw; not Mosflm MERIT 0-999).\n";
    if (lastgood > 1) {
        os << " SUGGLATT uses LASTGOOD rank " << lastgood
           << " (penalty jump before rank " << (lastgood + 1)
           << "), then highest symmetry; soft cut ≤ " << max_penalty
           << " is fallback.\n";
    } else {
        os << " SUGGLATT selects among PENALTY ≤ " << max_penalty
           << " (highest symmetry; no penalty jump detected).\n";
    }
    os << " No  PENALTY LATT     a        b        c"
       << "     alpha  beta gamma  Possible spacegroups\n\n";

    // Mosflm: DO I = NLAUEG,1,-1 — worst first, best last.
    // Widen PENALTY vs F77 I4: raw DPS criteria can exceed 9999.
    const std::ostream::fmtflags saved = os.flags();
    const std::streamsize old_prec = os.precision();
    for (std::size_t rev = order.size(); rev > 0; --rev) {
        const int no = static_cast<int>(rev);  // 1 = best
        const auto& c = candidates[order[rev - 1]];
        const int penalty = static_cast<int>(std::lround(c.criterion));
        const UnitCell cell = cell_from_dps_om(c.orientation);
        const std::string_view possible =
            possible_spacegroups_for_lattice(c.system_id);

        os << ' ' << std::setw(3) << no << ' ' << std::setw(7) << penalty
           << "  " << std::left << std::setw(2) << c.system_id << std::right;
        os << std::fixed << std::setprecision(2) << std::setw(9) << cell.a
           << std::setw(9) << cell.b << std::setw(9) << cell.c;
        os << "  " << std::setprecision(1) << std::setw(6) << cell.alpha_deg
           << std::setw(6) << cell.beta_deg << std::setw(6) << cell.gamma_deg
           << "   " << possible;
        if (lastgood > 1 && no == lastgood) {
            os << "  <--- LASTGOOD";
        }
        os << '\n';
        if (c.system_id == "hP") {
            // Mosflm FORMAT(65X,A) continuation for hexagonal P6 list.
            os << std::string(65, ' ') << kHexagonalP6 << '\n';
        }
    }
    os.flags(saved);
    os.precision(old_prec);

    const std::string_view sg_label = spgmin_label_for_lattice(selected_system);
    os << "\n\n Suggested Solution: " << suggested_rank << "  " << sg_label
       << "\n";
    if (lastgood > 1) {
        os << "           LASTGOOD:  " << lastgood
           << " (penalty jump before rank " << (lastgood + 1) << ")\n";
    }
    os << "            penalty:  "
       << static_cast<int>(std::lround(suggested_penalty)) << "\n"
       << "               cell:  " << std::fixed << std::setprecision(3)
       << selected_cell.a << "  " << selected_cell.b << "  " << selected_cell.c
       << "  " << std::setprecision(2) << selected_cell.alpha_deg << "  "
       << selected_cell.beta_deg << "  " << selected_cell.gamma_deg << "\n"
       << "           Symmetry: " << selected_system << " (" << sg_label
       << " / #" << symmetry_number_for_lattice(selected_system)
       << ", character " << selected_character << ")\n";
    os << std::flush;
}

}  // namespace mosflm::indexing
