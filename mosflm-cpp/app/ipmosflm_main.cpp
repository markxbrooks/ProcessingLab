/**
 * Thin ``ipmosflm`` CCP4-compatible name (Phase D / E).
 *
 * Default ProcessingLab release binary: CONTROL → C++ science path.
 * Does **not** link the Fortran amalgamation. Oracle:
 * ``make fortran-oracle`` (see ``../portingE.md``).
 */

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/diffraction/ProcessPipeline.hpp"

#if defined(MOSFLM_HAS_INDEXING)
#include "mosflm/indexing/ControlAutoIndex.hpp"
#endif

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
    // Fortran-compatible lead line + C++ core note (Phase E / E0).
    std::cout
        << "\n"
        << " ************ Mosflm version 7.4.0 for IP, CCD and PAD data"
           " (C++ core) ***********\n"
        << " ProcessingLab thin ipmosflm — CONTROL → C++ PROCESS/AUTOINDEX\n"
        << " Andrew Leslie and Harry Powell, MRC Laboratory of Molecular Biology\n"
        << " Oracle amalgamation: make fortran-oracle\n"
        << "--------------------------------------------------------------------------------\n";
}

void print_no_run_and_end() {
    std::cout << "\n\n *** WARNING **   NO \"RUN\" CARD GIVEN\n\n\n";
    std::cout << " *********** END OF PROCESSING *****************\n";
}

void print_end() {
    std::cout << "\n *********** END OF PROCESSING *****************\n";
}

}  // namespace

int main(int argc, char** argv) {
#if defined(MOSFLM_HAS_INDEXING)
    mosflm::indexing::register_control_autoindex();
#endif

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
            << "  Reads CONTROL from script or stdin (GO/RUN → C++ science).\n"
            << "  MOSFLM_IMAGE / CCP4 used when IMAGE/TEMPLATE unset.\n";
        return 0;
    }

    print_banner();

    try {
        mosflm::control::ControlSession session;
        if (hklout) {
            session.parse_line(std::string("HKLOUT ") + hklout->string());
        }

        std::optional<std::filesystem::path> img_override = image;
        if (!img_override) {
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
            if (r.ran_autoindex) {
                std::cout << "  AUTOINDEX: " << r.message << '\n';
                continue;
            }
            if (r.ran_strategy) {
                std::cout << "  STRATEGY: " << r.message << '\n';
                continue;
            }
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
            if (!session.state().extra_lattices.empty()) {
                std::cout << "  MULTILATTICE: "
                          << session.state().extra_lattices.size()
                          << " extra lattice(s) stored\n";
            }
        }

        if (session.needs_no_run_warning()) {
            print_no_run_and_end();
        } else {
            print_end();
        }
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "ipmosflm: " << ex.what() << '\n';
        return 1;
    }
}
