#pragma once

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/io/ExperimentModel.hpp"
#include "mosflm/io/Json.hpp"
#include "mosflm/io/RunMode.hpp"

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace mosflm::io {

/// Schema version of ``input.json`` and ``process.json``.
inline constexpr int kRunRecordSchemaVersion = 1;

/// How a run was invoked, for ``input.json``.
struct RunContext {
    std::string mode;  ///< ``"run"`` or ``"script"``
    RunMode run_mode{RunMode::PROCESSINGLAB};  ///< ``--run-mode``: parameter source
    std::vector<std::string> command_line;
    std::string control_text;  ///< CONTROL actually executed
    std::optional<std::filesystem::path> xds_inp;
    std::optional<std::filesystem::path> xinfo;
    std::optional<std::filesystem::path> xparm;
    std::optional<std::string> sweep;
};

/// Stage name for the step a ``GO`` will run (``"findspots"``,
/// ``"autoindex"``, ``"process"``, ``"strategy"``), from the pending flags.
[[nodiscard]] std::string pending_stage(const control::ControlState& state);

/// Stages a CONTROL text will run, in order, by parsing it without running
/// (``XPARMOUT`` is skipped: it writes a file when parsed).
[[nodiscard]] std::vector<std::string> planned_stages(const std::string& control_text);

/// Effective settings (indexing, integration, refinement) from a parsed
/// CONTROL state.
[[nodiscard]] Json settings_json(const control::ControlState& state);

/// ``input.json``: what ProcessingLab was asked to process, as resolved
/// before processing — images, imports, the experiment with per-field
/// provenance and overrides, supplied constraints and effective settings.
[[nodiscard]] Json input_record(
    const ExperimentModel& model,
    const control::ControlState& settings,
    const RunContext& context
);

/// ``process.json`` with lifecycle semantics, plus ``process.md``.
///
/// Created ``running`` with every planned stage ``pending``; each stage goes
/// ``running`` → ``completed`` (or ``failed``) with its results; the run ends
/// ``completed`` or ``failed`` with the error. Every change rewrites the
/// files atomically, so a crash leaves ``status: running`` and the stage
/// that was running. ``process.md`` is rendered from ``input.json`` and
/// ``process.json``, never from the log.
class ProcessRecord {
public:
    /// Writes ``input.json``, then ``process.json`` and ``process.md``.
    ProcessRecord(std::filesystem::path dir, Json input, const std::vector<std::string>& stages);

    /// A ``GO`` is starting (step still pending in ``state``).
    void begin_stage(const control::ControlState& state);
    /// It finished: record results and outputs.
    void end_stage(const control::ControlGoResult& result, const control::ControlState& state);
    /// Processing failed (the running stage, if any, fails too).
    void fail(const std::string& error);
    /// All stages done.
    void complete();
    /// A warning for ``process.json`` / ``process.md``.
    void add_warning(const std::string& warning);
    /// An output artifact (listed once, and only if it exists).
    void add_output(const std::filesystem::path& path, const std::string& kind);

    [[nodiscard]] const Json& document() const { return doc_; }

private:
    void write();
    Json& stage(const std::string& name);

    std::filesystem::path dir_;
    Json input_;
    Json doc_;
    std::string running_;
};

/// ``process.md`` from the two records; outputs inside ``dir`` are shown by
/// file name.
[[nodiscard]] std::string render_markdown(
    const Json& input, const Json& process, const std::filesystem::path& dir = {}
);

}  // namespace mosflm::io
