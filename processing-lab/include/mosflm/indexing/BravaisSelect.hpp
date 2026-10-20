#pragma once

#include "mosflm/geometry/UnitCell.hpp"
#include "mosflm/indexing/BravaisClassifier.hpp"

#include <cstddef>
#include <iosfwd>
#include <string>
#include <string_view>
#include <vector>

namespace mosflm::indexing {

/// Crystal-class rank for SUGGLATT-style selection (higher = more symmetry).
///
/// Matches Mosflm ICRYST ordering idea: cubic > hexagonal/rhombohedral >
/// tetragonal > orthorhombic > monoclinic > triclinic.
[[nodiscard]] int crystal_class_rank(std::string_view system_id);

/// Mosflm ``SPGMIN`` / DiffractionLab default space-group number for a
/// Bravais lattice code (e.g. ``tP`` → 75 = P4). Returns 1 (P1) if unknown.
[[nodiscard]] int symmetry_number_for_lattice(std::string_view system_id);

/// Mosflm ``SPLIST`` short label for the default space group of a lattice
/// (e.g. ``tP`` → ``P4``, ``cI`` → ``I23``).
[[nodiscard]] std::string_view spgmin_label_for_lattice(std::string_view system_id);

/// Mosflm ``SPLIST`` "Possible spacegroups" string for a Bravais lattice.
[[nodiscard]] std::string_view possible_spacegroups_for_lattice(
    std::string_view system_id
);

/// Mosflm DPSCHOOSE ``LASTGOOD``: 1-based rank (sorted ascending penalty)
/// of the last candidate before the largest relative MERIT jump.
///
/// Returns 1 when fewer than two candidates or no jump exceeds the Mosflm
/// baseline (0.1). MERIT is ``round(crit * 999 / max_crit)`` for jump math
/// only; printed PENALTY remains the raw DPS goof.
[[nodiscard]] int find_lastgood_rank(
    const std::vector<BravaisCandidate>& candidates
);

/// Select a Bravais candidate (Mosflm SUGGLATT / DPSCHOOSE penalty-jump).
///
/// 1. ``prefer_character`` (1–44) wins if in range.
/// 2. Else lowest-criterion match to ``prefer_system`` (e.g. ``"cI"``).
/// 3. Else Mosflm ``LASTGOOD`` jump window: among ranks ``1..LASTGOOD``
///    (when ``LASTGOOD > 1``), maximize crystal-class rank; tie-break lower
///    criterion then lower character.
/// 4. Else among candidates with ``criterion ≤ max_penalty``, maximize
///    crystal-class rank (soft-cut fallback).
/// 5. If none are under the penalty cut, fall back to global min criterion.
///
/// \return Index into ``candidates``, or ``npos`` if empty.
[[nodiscard]] std::size_t select_suggested_bravais(
    const std::vector<BravaisCandidate>& candidates,
    float max_penalty = 50.0f,
    int prefer_character = 0,
    std::string_view prefer_system = {}
);

/// Kabsch-style lattice metric check (≤3% axes, ≤3° angles vs ideal).
///
/// Used after IDXREF so SUGGLATT cannot promote mC/C2 on a skewed triclinic
/// cell. Triclinic ``aP`` always passes.
[[nodiscard]] bool lattice_metric_acceptable(
    std::string_view system_id,
    const UnitCell& cell,
    double max_axis_frac = 0.03,
    double max_angle_deg = 3.0
);

/// Post-IDXREF Bravais pick: keep refined aP unless a higher-symmetry
/// character passes :func:`lattice_metric_acceptable` and has a low criterion.
/// Prefers t/o/h/c over m when several pass. Honours ``prefer_*`` like
/// :func:`select_suggested_bravais`.
///
/// \param refined_cell Cell used for the metric gate (IDXREF result).
/// \return Index into ``candidates``, or ``npos`` if empty.
[[nodiscard]] std::size_t select_bravais_after_idxref(
    const std::vector<BravaisCandidate>& candidates,
    const UnitCell& refined_cell,
    float max_penalty = 50.0f,
    int prefer_character = 0,
    std::string_view prefer_system = {}
);

/// Write Mosflm-style lattice-character / penalty table (DPSCHOOSE format
/// without DI/SDCELL/FRACN — ``No PENALTY LATT a b c alpha beta gamma``).
///
/// Rows are sorted by ascending penalty then printed worst→best (``No`` =
/// 1-based rank after sort, matching Mosflm ``IRESORDER``). Marks the
/// SUGGLATT / LASTGOOD pick in the footer.
void write_bravais_penalty_table(
    std::ostream& os,
    const std::vector<BravaisCandidate>& candidates,
    int selected_character,
    std::string_view selected_system,
    const UnitCell& selected_cell,
    float max_penalty = 50.0f
);

}  // namespace mosflm::indexing
