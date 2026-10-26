#pragma once

#include <array>
#include <optional>
#include <string>
#include <string_view>

namespace mosflm::io {

/// A processing step of a run, in the order they run.
///
/// - ``FIND_SPOTS``: CONTROL ``FINDSPOTS`` (strong spots → ``.spt``).
/// - ``INDEXING``: ``AUTOINDEX`` (lattice, cell and orientation).
/// - ``DEFPIX``: the pixels integration ignores (untrusted regions,
///   excluded resolution, shaded / defective / hot pixels); runs at the
///   start of each ``PROCESS`` unless ``DEFPIX OFF``.
/// - ``INTEGRATE``: ``PROCESS`` (refine, predict, integrate → MTZ).
/// - ``CORRECT``: symmetry analysis of the integrated observations.
/// - ``STRATEGY``: ``STRATEGY`` (completeness estimate).
///
/// Records name a step by :func:`run_step_id` (``"defpix"``); logs and
/// messages by :func:`run_step_name` (``"RunStep.DEFPIX"``). DiffractionLab
/// mirrors this enum (``diffraction_lab.run_step.RunStep``).
enum class RunStep {
    FIND_SPOTS,
    INDEXING,
    DEFPIX,
    INTEGRATE,
    CORRECT,
    STRATEGY,
};

/// Every step, in run order.
inline constexpr std::array<RunStep, 6> kRunSteps{
    RunStep::FIND_SPOTS, RunStep::INDEXING, RunStep::DEFPIX, RunStep::INTEGRATE,
    RunStep::CORRECT, RunStep::STRATEGY};

/// ``"RunStep.FIND_SPOTS"``, … for logs and messages.
[[nodiscard]] inline const char* run_step_name(RunStep step) {
    switch (step) {
        case RunStep::FIND_SPOTS: return "RunStep.FIND_SPOTS";
        case RunStep::INDEXING: return "RunStep.INDEXING";
        case RunStep::DEFPIX: return "RunStep.DEFPIX";
        case RunStep::INTEGRATE: return "RunStep.INTEGRATE";
        case RunStep::CORRECT: return "RunStep.CORRECT";
        case RunStep::STRATEGY: break;
    }
    return "RunStep.STRATEGY";
}

/// The step's name in records (``process.json`` ``stages[].name``):
/// ``"find_spots"``, ``"indexing"``, ``"defpix"``, ``"integrate"``,
/// ``"strategy"``.
[[nodiscard]] inline const char* run_step_id(RunStep step) {
    switch (step) {
        case RunStep::FIND_SPOTS: return "find_spots";
        case RunStep::INDEXING: return "indexing";
        case RunStep::DEFPIX: return "defpix";
        case RunStep::INTEGRATE: return "integrate";
        case RunStep::CORRECT: return "correct";
        case RunStep::STRATEGY: break;
    }
    return "strategy";
}

/// The CONTROL keyword that makes the step pending (``DEFPIX`` runs inside
/// ``PROCESS``).
[[nodiscard]] inline const char* run_step_keyword(RunStep step) {
    switch (step) {
        case RunStep::FIND_SPOTS: return "FINDSPOTS";
        case RunStep::INDEXING: return "AUTOINDEX";
        case RunStep::DEFPIX:
        case RunStep::INTEGRATE: return "PROCESS";
        case RunStep::CORRECT: return "CORRECT";
        case RunStep::STRATEGY: break;
    }
    return "STRATEGY";
}

/// Parse a step: its id (``"defpix"``), name (``"RunStep.DEFPIX"`` or
/// ``"DEFPIX"``), any case, or a schema 1–2 stage name (``"findspots"``,
/// ``"autoindex"``, ``"process"``).
[[nodiscard]] inline std::optional<RunStep> parse_run_step(std::string_view s) {
    std::string v(s);
    for (auto& c : v) {
        c = static_cast<char>(c >= 'A' && c <= 'Z' ? c - 'A' + 'a' : c);
    }
    if (v.rfind("runstep.", 0) == 0) {
        v.erase(0, 8);
    }
    for (const RunStep step : kRunSteps) {
        if (v == run_step_id(step)) {
            return step;
        }
    }
    if (v == "findspots") return RunStep::FIND_SPOTS;
    if (v == "autoindex") return RunStep::INDEXING;
    if (v == "process") return RunStep::INTEGRATE;
    return std::nullopt;
}

}  // namespace mosflm::io
