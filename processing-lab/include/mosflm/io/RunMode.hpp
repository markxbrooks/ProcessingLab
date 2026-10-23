#pragma once

#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace mosflm::io {

/// Where a run's parameters come from (CLI ``--run-mode``).
///
/// - ``PROCESSINGLAB`` (default for ``ProcessingLab run``): image headers,
///   overridden per field by any ``--xds-inp`` / ``--xinfo`` imports.
/// - ``XDS``: an ``XDS.INP`` (``--xds-inp``, required); headers fill gaps.
/// - ``XIA2``: a Xia2 ``.xinfo`` (``--xinfo``, required); headers fill gaps.
/// - ``MOSFLM``: a MOSFLM-style CONTROL script (script mode, the default
///   there).
enum class RunMode {
    PROCESSINGLAB,
    XDS,
    MOSFLM,
    XIA2,
};

/// ``"RunMode.XDS"``, … for logs and records.
[[nodiscard]] inline const char* run_mode_name(RunMode mode) {
    switch (mode) {
        case RunMode::XDS: return "RunMode.XDS";
        case RunMode::MOSFLM: return "RunMode.MOSFLM";
        case RunMode::XIA2: return "RunMode.XIA2";
        case RunMode::PROCESSINGLAB: break;
    }
    return "RunMode.PROCESSINGLAB";
}

/// The ``--run-mode`` value: ``"processinglab"``, ``"xds"``, ``"mosflm"``,
/// ``"xia2"``.
[[nodiscard]] inline const char* run_mode_id(RunMode mode) {
    switch (mode) {
        case RunMode::XDS: return "xds";
        case RunMode::MOSFLM: return "mosflm";
        case RunMode::XIA2: return "xia2";
        case RunMode::PROCESSINGLAB: break;
    }
    return "processinglab";
}

/// Where this mode takes its parameters from, in words.
[[nodiscard]] inline const char* run_mode_parameter_source(RunMode mode) {
    switch (mode) {
        case RunMode::XDS: return "XDS.INP (--xds-inp); image headers fill unset fields";
        case RunMode::MOSFLM: return "MOSFLM CONTROL keywords (script); image headers fill unset fields";
        case RunMode::XIA2: return "Xia2 .xinfo (--xinfo); image headers fill unset fields";
        case RunMode::PROCESSINGLAB: break;
    }
    return "image headers, overridden per field by any --xds-inp / --xinfo";
}

/// Parse a ``--run-mode`` value (case-insensitive).
[[nodiscard]] inline std::optional<RunMode> parse_run_mode(std::string_view s) {
    std::string v(s);
    for (auto& c : v) {
        c = static_cast<char>(c >= 'A' && c <= 'Z' ? c - 'A' + 'a' : c);
    }
    if (v == "processinglab") return RunMode::PROCESSINGLAB;
    if (v == "xds") return RunMode::XDS;
    if (v == "mosflm") return RunMode::MOSFLM;
    if (v == "xia2") return RunMode::XIA2;
    return std::nullopt;
}

/// Check the inputs match the mode, so the parameter source is unambiguous.
///
/// :param mode: The run mode.
/// :param script_mode: True for a CONTROL script run (no ``run`` subcommand).
/// :param has_xds_inp: ``--xds-inp`` given.
/// :param has_xinfo: ``--xinfo`` given.
/// :param has_sweep: ``--sweep`` given.
/// :raises std::invalid_argument: With what to change, if they conflict.
inline void validate_run_mode(
    RunMode mode, bool script_mode, bool has_xds_inp, bool has_xinfo, bool has_sweep
) {
    const std::string name = run_mode_name(mode);
    switch (mode) {
        case RunMode::MOSFLM:
            if (!script_mode) {
                throw std::invalid_argument(
                    name + " takes its parameters from a CONTROL script: "
                    "ProcessingLab --run-mode mosflm [--image IMG] script.in (no 'run')"
                );
            }
            return;
        case RunMode::XDS:
            if (script_mode) {
                throw std::invalid_argument(name + " is for 'ProcessingLab run'; a CONTROL script is " +
                                            run_mode_name(RunMode::MOSFLM));
            }
            if (!has_xds_inp) {
                throw std::invalid_argument(name + " needs --xds-inp XDS.INP");
            }
            if (has_xinfo || has_sweep) {
                throw std::invalid_argument(name + " takes parameters from XDS.INP only; drop --xinfo/--sweep "
                                            "or use --run-mode xia2");
            }
            return;
        case RunMode::XIA2:
            if (script_mode) {
                throw std::invalid_argument(name + " is for 'ProcessingLab run'; a CONTROL script is " +
                                            run_mode_name(RunMode::MOSFLM));
            }
            if (!has_xinfo) {
                throw std::invalid_argument(name + " needs --xinfo FILE.xinfo");
            }
            if (has_xds_inp) {
                throw std::invalid_argument(name + " takes parameters from the .xinfo only; drop --xds-inp "
                                            "or use --run-mode xds");
            }
            return;
        case RunMode::PROCESSINGLAB:
            if (script_mode) {
                throw std::invalid_argument(name + " is 'ProcessingLab run'; a CONTROL script is " +
                                            run_mode_name(RunMode::MOSFLM));
            }
            return;
    }
}

}  // namespace mosflm::io
