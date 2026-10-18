/**
 * Thin ``ipmosflm`` CCP4-compatible name (Phase D / D6).
 *
 * Default ProcessingLab release binary: CONTROL keyword subset → C++ PROCESS.
 * Does **not** link the Fortran amalgamation. Rebuild the Autotools oracle with
 * ``make MOSFLM_FORTRAN_ORACLE=1`` (see ``../portingD.md``).
 *
 * Usage::
 *
 *   ipmosflm [HKLOUT out.mtz] [--image img] [script.in]
 *   ipmosflm HKLOUT out.mtz < script.in
 *   ipmosflm          # banner + read CONTROL from stdin until EXIT
 */

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/diffraction/ProcessPipeline.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace {

std::filesystem::path default_insulin_img() {
    if (const char* env = std::getenv("MOSFLM_IMAGE")) {
        return env;
    }
    if (const char* ccp4 = std::getenv("CCP4")) {
        const auto p =
            std::filesystem::path(ccp4) / "examples/data/insulin_1_001.img";
        if (std::filesystem::exists(p)) {
            return p;
        }
    }
    return "/Applications/ccp4-9/examples/data/insulin_1_001.img";
}

void print_banner() {
    std::cout
        << " ProcessingLab / Mosflm C++ core 0.1.0\n"
        << "  (thin ipmosflm — CONTROL subset; no Fortran amalgamation)\n"
        << "  Full oracle: make MOSFLM_FORTRAN_ORACLE=1 -C mosflm ipmosflm\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::optional<std::filesystem::path> image;
    std::optional<std::filesystem::path> hklout;
    std::optional<std::filesystem::path> script;
    bool want_help = false;

    for (int i = 1; i < argc; ++i) {
        const std::string_view a = argv[i];
        if (a == "-h" || a == "--help" || a == "help") {
            want_help = true;
            break;
        }
        if (a == "--image") {
            if (i + 1 >= argc) {
                std::cerr << "ipmosflm: --image needs a path\n";
                return 2;
            }
            image = argv[++i];
            continue;
        }
        if (a == "HKLOUT" || a == "hklout") {
            if (i + 1 >= argc) {
                std::cerr << "ipmosflm: HKLOUT needs a path\n";
                return 2;
            }
            hklout = argv[++i];
            continue;
        }
        if (a == "version" || a == "--version") {
            print_banner();
            return 0;
        }
        // First non-option positional: script path.
        if (!script) {
            script = argv[i];
            continue;
        }
        std::cerr << "ipmosflm: unexpected argument '" << a << "'\n";
        return 2;
    }

    if (want_help) {
        std::cerr
            << "Usage: ipmosflm [HKLOUT out.mtz] [--image img] [script.in]\n"
            << "  Reads CONTROL from script or stdin (GO → C++ PROCESS).\n"
            << "  MOSFLM_IMAGE / CCP4 used when IMAGE/TEMPLATE unset.\n";
        return 0;
    }

    print_banner();

    try {
        mosflm::control::ControlSession session;
        if (hklout) {
            session.parse_line(
                std::string("HKLOUT ") + hklout->string()
            );
        }

        std::optional<std::filesystem::path> img_override = image;
        if (!img_override) {
            // Prefer env / CCP4 when scripts omit IMAGE (known-indexed corpus).
            const auto fallback = default_insulin_img();
            if (std::filesystem::exists(fallback)) {
                img_override = fallback;
            }
        }

        std::vector<mosflm::control::ControlGoResult> runs;
        if (script) {
            runs = session.run_script(*script, img_override);
        } else {
            runs = session.run_stream(std::cin, img_override);
        }

        const auto out_mtz = hklout ? *hklout : session.state().hklout;
        std::size_t n_proc = 0;
        for (const auto& r : runs) {
            if (!r.ran_process) {
                continue;
            }
            ++n_proc;
            mosflm::diffraction::ProcessPipeline::write_mtz(
                out_mtz, r.crystal, r.process
            );
            std::cout << "  PROCESS: image=" << r.image_used
                      << " preds=" << r.process.predictions.size()
                      << " integrated=" << r.process.mtz_reflections.size()
                      << " → " << out_mtz << '\n';
        }
        if (n_proc == 0 && !session.exit_requested()) {
            // EXIT-only smoke (no PROCESS) is success.
            if (runs.empty()) {
                std::cout << "  (no PROCESS GO)\n";
            }
        }
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "ipmosflm: " << ex.what() << '\n';
        return 1;
    }
}
