/**
 * ``IntegrationLab`` — ProcessingLab C++ CONTROL driver (Phase D / E).
 *
 * Default release binary: CONTROL → C++ science path.
 * Does **not** link the Fortran amalgamation. F77 oracle remains ``ipmosflm``:
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
#include <iomanip>
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
        << " IntegrationLab — CONTROL → C++ PROCESS/AUTOINDEX\n"
        << " Andrew Leslie and Harry Powell, MRC Laboratory of Molecular Biology\n"
        << " F77 oracle (ipmosflm): make fortran-oracle\n"
        << "--------------------------------------------------------------------------------\n";
}

void print_no_run_and_end() {
    std::cout << "\n\n *** WARNING **   NO \"RUN\" CARD GIVEN\n\n\n";
    std::cout << " *********** END OF PROCESSING *****************\n";
}

void print_end() {
    std::cout << "\n *********** END OF PROCESSING *****************\n";
}

/// Machine-readable AUTOINDEX dump for DiffractionLab (and CLI users).
void emit_autoindex_dump(const mosflm::control::ControlState& state) {
    if (!state.has_cell || !state.has_matrix) {
        return;
    }
    const auto& c = state.cell;
    const auto& A = state.amat;
    std::cout << std::setprecision(10);
    std::cout << "INTEGRATIONLAB_CELL " << c.a << " " << c.b << " " << c.c
              << " " << c.alpha_deg << " " << c.beta_deg << " " << c.gamma_deg
              << '\n';
    std::cout << "INTEGRATIONLAB_MATRIX";
    for (int r = 0; r < 3; ++r) {
        for (int col = 0; col < 3; ++col) {
            std::cout << " " << A(r, col);
        }
    }
    std::cout << '\n';
    if (!state.autoindex_system.empty()) {
        std::cout << "INTEGRATIONLAB_SYSTEM " << state.autoindex_system << '\n';
    }
    if (state.autoindex_character > 0) {
        std::cout << "INTEGRATIONLAB_CHARACTER " << state.autoindex_character
                  << '\n';
    }
    if (state.autoindex_n_spots > 0) {
        std::cout << "INTEGRATIONLAB_NSPOTS " << state.autoindex_n_spots << '\n';
    }

    // Sidecar next to cwd for tools that prefer files over stdout.
    try {
        std::ofstream out("integrationlab_index.txt");
        if (out) {
            out << std::setprecision(10);
            out << "CELL " << c.a << " " << c.b << " " << c.c << " "
                << c.alpha_deg << " " << c.beta_deg << " " << c.gamma_deg
                << '\n';
            out << "MATRIX";
            for (int r = 0; r < 3; ++r) {
                for (int col = 0; col < 3; ++col) {
                    out << " " << A(r, col);
                }
            }
            out << '\n';
            if (!state.autoindex_system.empty()) {
                out << "SYSTEM " << state.autoindex_system << '\n';
            }
            if (state.autoindex_character > 0) {
                out << "CHARACTER " << state.autoindex_character << '\n';
            }
            if (state.autoindex_n_spots > 0) {
                out << "NSPOTS " << state.autoindex_n_spots << '\n';
            }
        }
    } catch (...) {
        // Best-effort sidecar; stdout dump is enough for DiffractionLab.
    }
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
                std::cerr << "IntegrationLab: --image needs a path\n";
                return 2;
            }
            image = argv[++i];
            continue;
        }
        if (a == "HKLOUT" || a == "hklout") {
            if (i + 1 >= argc) {
                std::cerr << "IntegrationLab: HKLOUT needs a path\n";
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
        std::cerr << "IntegrationLab: unexpected argument '" << a << "'\n";
        return 2;
    }

    if (want_help) {
        std::cerr
            << "Usage: IntegrationLab [HKLOUT out.mtz] [--image img] [script.in]\n"
            << "  Reads CONTROL from script or stdin (GO/RUN → C++ science).\n"
            << "  MOSFLM_IMAGE / CCP4 used when IMAGE/TEMPLATE unset.\n"
            << "  F77 ipmosflm oracle: make fortran-oracle\n";
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
                emit_autoindex_dump(session.state());
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
        std::cerr << "IntegrationLab: " << ex.what() << '\n';
        return 1;
    }
}
