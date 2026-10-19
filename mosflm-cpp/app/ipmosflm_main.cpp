/**
 * ``ProcessingLab`` — ProcessingLab C++ CONTROL driver (Phase D / E).
 *
 * Default release binary: CONTROL → C++ science path.
 * Does **not** link the Fortran amalgamation. F77 oracle remains ``ipmosflm``:
 * ``make fortran-oracle`` (see ``../portingE.md``).
 */

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/diffraction/ProcessPipeline.hpp"
#include "mosflm/io/ImportedSweepGeometry.hpp"
#include "mosflm/io/XdsInpImporter.hpp"
#include "mosflm/io/XinfoImporter.hpp"
#include "mosflm/io/XparmIO.hpp"

#if defined(MOSFLM_HAS_INDEXING)
#include "mosflm/indexing/ControlAutoIndex.hpp"
#endif

#include <cmath>
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
        << "============================================================\n"
        << " ProcessingLab v0.1.0\n"
        << "------------------------------------------------------------\n"
        << " Diffraction Processing & Autoindexing\n"
        << " Data: CCD | CBF | HDF5\n"
        << "------------------------------------------------------------\n"
        << " Developed by Mark Brooks 2026\n"
        << "============================================================\n"
        << std::flush;
}

void print_no_run_and_end() {
    std::cout << "\n\n *** WARNING **   NO \"RUN\" CARD GIVEN\n\n\n";
    std::cout << " ============= END OF PROCESSING ================\n";
}

void print_end() {
    std::cout << " ============= END OF PROCESSING ================\n";
}

/// Machine-readable AUTOINDEX dump for DiffractionLab (and CLI users).
void emit_autoindex_dump(const mosflm::control::ControlState& state) {
    if (!state.has_cell || !state.has_matrix) {
        return;
    }
    const auto& c = state.cell;
    const auto& A = state.amat;
    std::cout << std::setprecision(10);
    std::cout << "PROCESSINGLAB_CELL " << c.a << " " << c.b << " " << c.c
              << " " << c.alpha_deg << " " << c.beta_deg << " " << c.gamma_deg
              << '\n';
    std::cout << "PROCESSINGLAB_MATRIX";
    for (int r = 0; r < 3; ++r) {
        for (int col = 0; col < 3; ++col) {
            std::cout << " " << A(r, col);
        }
    }
    std::cout << '\n';
    if (!state.autoindex_system.empty()) {
        std::cout << "PROCESSINGLAB_SYSTEM " << state.autoindex_system << '\n';
    }
    if (state.autoindex_character > 0) {
        std::cout << "PROCESSINGLAB_CHARACTER " << state.autoindex_character
                  << '\n';
    }
    if (state.autoindex_n_spots > 0) {
        std::cout << "PROCESSINGLAB_NSPOTS " << state.autoindex_n_spots << '\n';
    }
    if (state.autoindex_indexed_frac > 0.0 || state.autoindex_n_spots > 0) {
        std::cout << "PROCESSINGLAB_INDEXED "
                  << static_cast<int>(
                         std::lround(100.0 * state.autoindex_indexed_frac)
                     )
                  << '\n';
        std::cout << std::setprecision(4);
        std::cout << "PROCESSINGLAB_RMSD_HKL " << state.autoindex_rmsd_hkl
                  << '\n';
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
            out << "INDEXED "
                << static_cast<int>(
                       std::lround(100.0 * state.autoindex_indexed_frac)
                   )
                << '\n';
            out << std::setprecision(4);
            out << "RMSD_HKL " << state.autoindex_rmsd_hkl << '\n';
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
    std::optional<std::filesystem::path> xinfo_path;
    std::optional<std::filesystem::path> xds_inp_path;
    std::optional<std::filesystem::path> xparm_path;
    bool want_help = false;

    for (int i = 1; i < argc; ++i) {
        const std::string_view a = argv[i];
        if (a == "-h" || a == "--help" || a == "help") {
            want_help = true;
            break;
        }
        if (a == "--image") {
            if (i + 1 >= argc) {
                std::cerr << "ProcessingLab: --image needs a path\n";
                return 2;
            }
            image = argv[++i];
            continue;
        }
        if (a == "--xinfo") {
            if (i + 1 >= argc) {
                std::cerr << "ProcessingLab: --xinfo needs a path\n";
                return 2;
            }
            xinfo_path = argv[++i];
            continue;
        }
        if (a == "--xds-inp") {
            if (i + 1 >= argc) {
                std::cerr << "ProcessingLab: --xds-inp needs a path\n";
                return 2;
            }
            xds_inp_path = argv[++i];
            continue;
        }
        if (a == "--xparm" || a == "--gxparm") {
            if (i + 1 >= argc) {
                std::cerr << "ProcessingLab: --xparm needs a path\n";
                return 2;
            }
            xparm_path = argv[++i];
            continue;
        }
        if (a == "HKLOUT" || a == "hklout") {
            if (i + 1 >= argc) {
                std::cerr << "ProcessingLab: HKLOUT needs a path\n";
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
        std::cerr << "ProcessingLab: unexpected argument '" << a << "'\n";
        return 2;
    }

    if (want_help) {
        std::cerr
            << "Usage: ProcessingLab [HKLOUT out.mtz] [--image img]\n"
            << "         [--xinfo file.xinfo] [--xds-inp XDS.INP]\n"
            << "         [--xparm XPARM.XDS] [script.in]\n"
            << "  Reads CONTROL from script or stdin (GO/RUN → C++ science).\n"
            << "  --xds-inp then --xinfo then --xparm fill CONTROL before script.\n"
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

        // XDS.INP first, then .xinfo (xinfo wins for paths / range / λ).
        if (xds_inp_path) {
            mosflm::io::apply_imported_geometry(
                session, mosflm::io::XdsInpImporter::parse_file(*xds_inp_path)
            );
            std::cout << "  Imported XDS.INP: " << *xds_inp_path << '\n';
        }
        if (xinfo_path) {
            mosflm::io::apply_imported_geometry(
                session, mosflm::io::XinfoImporter::parse_file(*xinfo_path)
            );
            std::cout << "  Imported xinfo: " << *xinfo_path << '\n';
        }
        if (xparm_path) {
            mosflm::io::apply_xparm(
                session, mosflm::io::read_xparm(*xparm_path)
            );
            std::cout << "  Imported XPARM: " << *xparm_path << '\n';
        }

        std::optional<std::filesystem::path> img_override = image;
        const bool imported_images = xinfo_path || xds_inp_path ||
            !session.state().image_path.empty() ||
            !session.state().image_template.empty();
        if (!img_override && !imported_images) {
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
            if (r.ran_findspots) {
                std::cout << "  FINDSPOTS: " << r.message << '\n';
                continue;
            }
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
                out_mtz,
                r.crystal,
                r.process,
                session.state().wavelength_A
            );
            std::cout << "  PROCESS: image=" << r.image_used;
            if (r.process_last > r.process_first) {
                std::cout << " frames=" << r.process_first << "-"
                          << r.process_last;
            } else if (r.process_first > 0) {
                std::cout << " frame=" << r.process_first;
            }
            std::cout << " preds=" << r.process.predictions.size()
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
        std::cerr << "ProcessingLab: " << ex.what() << '\n';
        return 1;
    }
}
