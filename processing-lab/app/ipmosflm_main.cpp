/**
 * ``ProcessingLab`` — ProcessingLab C++ CONTROL driver (Phase D / E).
 *
 * Default release binary: CONTROL → C++ science path.
 * Does **not** link the Fortran amalgamation. F77 oracle remains ``ipmosflm``:
 * ``make fortran-oracle`` (see ``../portingE.md``).
 */

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/diffraction/ProcessPipeline.hpp"
#include "mosflm/io/ExperimentModel.hpp"
#include "mosflm/io/ImportedSweepGeometry.hpp"
#include "mosflm/io/RunMode.hpp"
#include "mosflm/io/RunRecord.hpp"
#include "mosflm/io/XdsInpImporter.hpp"
#include "mosflm/io/XinfoImporter.hpp"
#include "mosflm/io/XparmIO.hpp"

#if defined(MOSFLM_HAS_INDEXING)
#include "mosflm/indexing/ControlAutoIndex.hpp"
#endif

#include <algorithm>
#include <cctype>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

namespace {

/// Format the CONTROL input as a block at the top of the log.
std::string control_input_block(const std::string& text) {
    std::string out = "\n ***** CONTROL INPUT *****\n\n" + text;
    if (!text.empty() && text.back() != '\n') {
        out += '\n';
    }
    out += "\n ***** END OF CONTROL INPUT *****\n\n";
    return out;
}

/// Copies everything written to stdout / stderr (C++ streams, C ``printf``
/// from the DPS indexer, …) into a log file (default ``process.log``) while
/// still showing it on the terminal.
///
/// Descriptors 1 and 2 are redirected into pipes; one thread per pipe copies
/// to the original descriptor and to the log. ``close()`` restores them and
/// drains the pipes; it also runs from ``atexit`` so output printed before an
/// ``exit()`` reaches the terminal.
class ProcessLog {
public:
    ProcessLog(const std::filesystem::path& path, const std::string& header) {
        log_fd_ = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (log_fd_ < 0) {
            throw std::runtime_error(
                "cannot open log file '" + path.string() + "'"
            );
        }
        write_all(log_fd_, header.data(), header.size());
        std::cout.flush();
        std::cerr.flush();
        std::fflush(stdout);
        std::fflush(stderr);
        redirect(STDOUT_FILENO, out_);
        redirect(STDERR_FILENO, err_);
        // A pipe makes stdio fully buffered: keep stdout line-buffered so
        // prompts and progress still appear promptly.
        std::setvbuf(stdout, nullptr, _IOLBF, 0);
        std::setvbuf(stderr, nullptr, _IONBF, 0);
        active_ = this;
        std::atexit([] {
            if (active_ != nullptr) {
                active_->close();
            }
        });
    }

    ProcessLog(const ProcessLog&) = delete;
    ProcessLog& operator=(const ProcessLog&) = delete;

    ~ProcessLog() { close(); }

    /// Restore stdout / stderr, drain the pipes and close the log.
    void close() {
        if (log_fd_ < 0) {
            return;
        }
        active_ = nullptr;
        std::cout.flush();
        std::cerr.flush();
        std::fflush(stdout);
        std::fflush(stderr);
        restore(STDOUT_FILENO, out_);
        restore(STDERR_FILENO, err_);
        ::close(log_fd_);
        log_fd_ = -1;
    }

private:
    struct Channel {
        int saved_fd{-1};
        int read_fd{-1};
        std::thread pump;
    };

    static void write_all(int fd, const char* p, std::size_t n) {
        while (n > 0) {
            const ssize_t w = ::write(fd, p, n);
            if (w <= 0) {
                return;
            }
            p += w;
            n -= static_cast<std::size_t>(w);
        }
    }

    void redirect(int fd, Channel& ch) {
        int fds[2];
        if (::pipe(fds) != 0) {
            throw std::runtime_error("cannot create log pipe");
        }
        ch.saved_fd = ::dup(fd);
        ch.read_fd = fds[0];
        ::dup2(fds[1], fd);
        ::close(fds[1]);
        ch.pump = std::thread([this, &ch] {
            char buf[8192];
            ssize_t n = 0;
            while ((n = ::read(ch.read_fd, buf, sizeof buf)) > 0) {
                write_all(ch.saved_fd, buf, static_cast<std::size_t>(n));
                const std::lock_guard<std::mutex> lock(log_mutex_);
                write_all(log_fd_, buf, static_cast<std::size_t>(n));
            }
        });
    }

    static void restore(int fd, Channel& ch) {
        if (ch.saved_fd < 0) {
            return;
        }
        ::dup2(ch.saved_fd, fd);  // last write end closed → pump sees EOF
        if (ch.pump.joinable()) {
            ch.pump.join();
        }
        ::close(ch.read_fd);
        ::close(ch.saved_fd);
        ch.saved_fd = -1;
    }

    static inline ProcessLog* active_ = nullptr;
    int log_fd_{-1};
    std::mutex log_mutex_;
    Channel out_;
    Channel err_;
};

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

    // Sidecar (output directory, else cwd) for tools that prefer files.
    try {
        std::ofstream out(
            state.output_directory.empty()
                ? std::filesystem::path("index.log")
                : state.output_directory / "index.log"
        );
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

/// Write MTZ for each PROCESS and report every GO (FINDSPOTS, AUTOINDEX, …).
void report_runs(
    const mosflm::control::ControlSession& session,
    const std::vector<mosflm::control::ControlGoResult>& runs,
    const std::filesystem::path& out_mtz
) {
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
}

/// Parse CONTROL into ``session`` without running any ``GO``, to record the
/// effective settings before processing. ``XPARMOUT`` is skipped (it writes a
/// file when parsed); parsing stops at ``EXIT``.
void dry_parse(mosflm::control::ControlSession& session, const std::string& control) {
    std::istringstream in(control);
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream toks(line.substr(0, line.find('!')));
        std::string key;
        toks >> key;
        for (auto& c : key) {
            c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        }
        if (key.empty() || key == "GO" || key == "RUN" || key == "XPARMOUT") {
            continue;
        }
        if (key == "EXIT" || key == "QUIT" || key == "END") {
            break;
        }
        session.parse_line(line);
    }
}

/// Wire ``record`` to every GO of ``session``.
void observe(mosflm::control::ControlSession& session, mosflm::io::ProcessRecord& record) {
    session.set_go_observer({
        [&record](const mosflm::control::ControlState& st) { record.begin_stage(st); },
        [&record](const mosflm::control::ControlGoResult& r, const mosflm::control::ControlState& st) {
            record.end_stage(r, st);
        },
    });
}

void print_run_usage() {
    std::cout
        << "Usage: ProcessingLab run --image PATH [--image PATH ...]\n"
        << "         [--range FIRST:LAST] [--output-dir DIR]\n"
        << "         [--xds-inp XDS.INP] [--xinfo FILE.xinfo [--sweep NAME]]\n"
        << "         [--run-mode processinglab|xds|xia2] [--help] [--version]\n"
        << "\n"
        << "  Finds spots, autoindexes and integrates each sweep. Metadata comes\n"
        << "  from the image headers; XDS.INP, then .xinfo, then the command line\n"
        << "  override it field by field, and the log shows each value's source.\n"
        << "  --image      HDF5/NeXus master, or one frame of a CBF/SMV sweep\n"
        << "               (repeatable: one sweep per image, processed in turn).\n"
        << "               Optional when --xinfo or --xds-inp names the images.\n"
        << "  --range      Inclusive frame range (default: every frame).\n"
        << "  --output-dir New directory for all outputs (must not exist). Default:\n"
        << "               processing_YYYYMMDD_HHMMSS in the current directory.\n"
        << "               Holds input.json (resolved inputs and provenance),\n"
        << "               process.json (stage status and results), process.md,\n"
        << "               process.log, process.mtz, XPARM.XDS and GXPARM.XDS.\n"
        << "  --sweep      Which .xinfo sweep, if it has several.\n"
        << "  --run-mode   Where the parameters come from (recorded in the log and\n"
        << "               input.json):\n"
        << "                 processinglab  image headers, overridden by any imports\n"
        << "                                (RunMode.PROCESSINGLAB, the default)\n"
        << "                 xds            XDS.INP; needs --xds-inp, no --xinfo\n"
        << "                                (RunMode.XDS)\n"
        << "                 xia2           .xinfo; needs --xinfo, no --xds-inp\n"
        << "                                (RunMode.XIA2)\n"
        << "               mosflm (RunMode.MOSFLM) is a CONTROL script run:\n"
        << "               ProcessingLab --run-mode mosflm [--image IMG] script.in\n";
}

/// ``ProcessingLab run``: image-first processing into a fresh directory.
int run_command(int argc, char** argv) {
    namespace fs = std::filesystem;
    std::vector<std::string> command_line{"ProcessingLab", "run"};
    for (int i = 0; i < argc; ++i) {
        command_line.emplace_back(argv[i]);
    }
    mosflm::io::ExperimentInputs base;
    std::vector<fs::path> images;
    std::optional<fs::path> output_dir;
    mosflm::io::RunMode run_mode = mosflm::io::RunMode::PROCESSINGLAB;
    const auto need = [&](int& i, std::string_view opt) -> const char* {
        if (i + 1 >= argc) {
            throw std::invalid_argument(std::string(opt) + " needs a value");
        }
        return argv[++i];
    };
    try {
        for (int i = 0; i < argc; ++i) {
            const std::string_view a = argv[i];
            if (a == "-h" || a == "--help") {
                print_run_usage();
                return 0;
            }
            if (a == "--version") {
                print_banner();
                return 0;
            }
            if (a == "--image") {
                images.emplace_back(need(i, a));
            } else if (a == "--range") {
                const std::string v = need(i, a);
                const auto colon = v.find(':');
                std::size_t used = 0;
                int first = 0;
                int last = 0;
                try {
                    first = std::stoi(v.substr(0, colon), &used);
                    if (used != colon) {
                        throw std::invalid_argument(v);
                    }
                    const std::string tail = v.substr(colon + 1);
                    last = std::stoi(tail, &used);
                    if (colon == std::string::npos || used != tail.size()) {
                        throw std::invalid_argument(v);
                    }
                } catch (const std::exception&) {
                    throw std::invalid_argument("--range needs FIRST:LAST, got '" + v + "'");
                }
                if (first < 1 || last < first) {
                    throw std::invalid_argument("--range needs 1 <= FIRST <= LAST, got '" + v + "'");
                }
                base.range = std::pair{first, last};
            } else if (a == "--output-dir") {
                output_dir = need(i, a);
            } else if (a == "--xds-inp") {
                base.xds_inp = fs::absolute(need(i, a));
            } else if (a == "--xinfo") {
                base.xinfo = fs::absolute(need(i, a));
            } else if (a == "--sweep") {
                base.sweep = need(i, a);
            } else if (a == "--run-mode") {
                const std::string v = need(i, a);
                const auto rm = mosflm::io::parse_run_mode(v);
                if (!rm) {
                    throw std::invalid_argument("--run-mode must be processinglab, xds, mosflm or xia2, got '" + v + "'");
                }
                run_mode = *rm;
            } else {
                throw std::invalid_argument("unexpected argument '" + std::string(a) + "' (see run --help)");
            }
        }
        mosflm::io::validate_run_mode(run_mode, /*script_mode=*/false, base.xds_inp.has_value(),
                                      base.xinfo.has_value(), base.sweep.has_value());
    } catch (const std::exception& ex) {
        std::cerr << "ProcessingLab run: " << ex.what() << '\n';
        return 2;
    }

    // Resolve every sweep before creating anything, so a bad input fails
    // without leaving an empty output directory.
    std::vector<mosflm::io::ExperimentModel> models;
    try {
        if (images.empty()) {
            models.push_back(mosflm::io::resolve_experiment(base));
        }
        for (const auto& img : images) {
            auto in = base;
            in.image = fs::absolute(img);
            models.push_back(mosflm::io::resolve_experiment(in));
        }
    } catch (const std::exception& ex) {
        std::cerr << "ProcessingLab run: " << ex.what() << '\n';
        return 2;
    }

    fs::path out_root;
    try {
        out_root = mosflm::io::create_output_directory(output_dir, fs::current_path(), std::time(nullptr));
    } catch (const std::exception& ex) {
        std::cerr << "ProcessingLab run: " << ex.what() << '\n';
        return 2;
    }

    int status = 0;
    for (std::size_t k = 0; k < models.size(); ++k) {
        const auto& m = models[k];
        // One sweep per directory: the run directory itself, or sweep_N in it.
        fs::path dir = out_root;
        if (models.size() > 1) {
            dir = out_root / ("sweep_" + std::to_string(k + 1));
            fs::create_directory(dir);
        }
        const auto [f0, f1] = m.frame_range.value;
        std::ostringstream control;
        for (const auto& line : m.control_lines()) {
            control << line << '\n';
        }
        if (m.cell.set()) {
            // DMAX also bounds the autoindexing search: the 200 Å default lost
            // 76c1b's 187 Å axis, while 320 Å skewed thaumatin's cell. With a
            // known cell, allow 1.7× its longest edge.
            const double edge = std::max({m.cell.value[0], m.cell.value[1], m.cell.value[2]});
            control << "DMAX " << std::max(200.0, std::ceil(1.7 * edge)) << '\n';
        }
        control << "HKLOUT " << (dir / "process.mtz").string() << '\n'
                << "FINDSPOTS THRESHOLD 5.0 MAXSPOTS 1500 FIND " << f0
                << " FILE spots.spt MULTIFRAME RANGE " << f0 << ' ' << std::min(f0 + 49, f1) << '\n'
                << "GO\n"
                << "AUTOINDEX DPS THRESHOLD 15\n"
                << "GO\n"
                << "PROCESS " << f0 << " TO " << f1;
        if (m.oscillation_deg.set() && m.oscillation_deg.source != mosflm::io::MetadataSource::ImageHeader) {
            control << " ANGLE " << m.oscillation_deg.value;
        }
        control << "\nGO\nEXIT\n";

        std::ostringstream header;
        header << "\n Run mode:        " << mosflm::io::run_mode_name(run_mode) << "  (parameters from "
               << mosflm::io::run_mode_parameter_source(run_mode) << ")\n";
        m.write_summary(header);
        header << "  Output:        " << dir.string() << '\n';

        std::optional<mosflm::io::ProcessRecord> record;
        try {
            ProcessLog plog(dir / "process.log", control_input_block(control.str()));
            print_banner();
            std::cout << header.str();
            // input.json: the effective inputs, settled before processing.
            mosflm::control::ControlSession dry;
            dry.set_output_directory(dir);
            dry_parse(dry, control.str());
            const mosflm::io::RunContext ctx{
                "run", run_mode, command_line, control.str(), base.xds_inp, base.xinfo, std::nullopt, base.sweep};
            record.emplace(dir, mosflm::io::input_record(m, dry.state(), ctx),
                           mosflm::io::planned_stages(control.str()));
            mosflm::control::ControlSession session;
            session.set_output_directory(dir);
            observe(session, *record);
            std::optional<fs::path> img_override;
            if (m.hdf5) {
                img_override = m.image.value;
            }
            std::istringstream in(control.str());
            const auto runs = session.run_stream(in, img_override, "");
            report_runs(session, runs, dir / "process.mtz");
            record->add_output(dir / "process.mtz", "mtz");
            record->add_output(dir / "index.log", "index_summary");
            record->add_output(dir / "process.log", "log");
            record->complete();
            print_end();
        } catch (const std::exception& ex) {
            std::cerr << "ProcessingLab run: " << ex.what() << '\n';
            if (record) {
                record->fail(ex.what());
            }
            status = 1;
        }
    }
    std::cout << "\n Output directory: " << out_root.string() << '\n';
    return status;
}

}  // namespace

int main(int argc, char** argv) {
#if defined(MOSFLM_HAS_INDEXING)
    mosflm::indexing::register_control_autoindex();
#endif
    if (argc > 1 && std::string_view(argv[1]) == "run") {
        return run_command(argc - 2, argv + 2);
    }

    std::optional<std::filesystem::path> image;
    std::optional<std::filesystem::path> hklout;
    std::optional<std::filesystem::path> script;
    std::optional<std::filesystem::path> xinfo_path;
    std::optional<std::filesystem::path> xds_inp_path;
    std::optional<std::filesystem::path> xparm_path;
    std::optional<std::filesystem::path> log_path = "process.log";
    std::optional<std::filesystem::path> output_dir;
    // CONTROL keywords are MOSFLM's: a script run is RunMode.MOSFLM.
    mosflm::io::RunMode run_mode = mosflm::io::RunMode::MOSFLM;
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
        if (a == "--log") {
            if (i + 1 >= argc) {
                std::cerr << "ProcessingLab: --log needs a path\n";
                return 2;
            }
            log_path = argv[++i];
            continue;
        }
        if (a == "--no-log") {
            log_path.reset();
            continue;
        }
        if (a == "--output-dir") {
            if (i + 1 >= argc) {
                std::cerr << "ProcessingLab: --output-dir needs a path\n";
                return 2;
            }
            output_dir = argv[++i];
            continue;
        }
        if (a == "--run-mode") {
            const auto rm = i + 1 < argc ? mosflm::io::parse_run_mode(argv[i + 1]) : std::nullopt;
            if (!rm) {
                std::cerr << "ProcessingLab: --run-mode must be processinglab, xds, mosflm or xia2\n";
                return 2;
            }
            ++i;
            try {
                mosflm::io::validate_run_mode(*rm, /*script_mode=*/true, false, false, false);
            } catch (const std::exception& ex) {
                std::cerr << "ProcessingLab: " << ex.what() << '\n';
                return 2;
            }
            run_mode = *rm;
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
            << "Usage: ProcessingLab run --image PATH ...   (see run --help)\n"
            << "       ProcessingLab [HKLOUT out.mtz] [--image img]\n"
            << "         [--xinfo file.xinfo] [--xds-inp XDS.INP]\n"
            << "         [--xparm XPARM.XDS] [--log process.log | --no-log]\n"
            << "         [--output-dir DIR] [--run-mode mosflm] [script.in]\n"
            << "  A CONTROL script run is RunMode.MOSFLM (parameters from MOSFLM\n"
            << "  keywords); --run-mode mosflm states it explicitly.\n"
            << "  --output-dir: new directory (must not exist) for the log, MTZ,\n"
            << "  spot list, XPARM.XDS, GXPARM.XDS and the run records input.json,\n"
            << "  process.json and process.md; script paths stay relative to the\n"
            << "  current directory.\n"
            << "  Output is also written to process.log (--log path to rename,\n"
            << "  --no-log to disable), including the CONTROL input.\n"
            << "  Reads CONTROL from script or stdin (GO/RUN → C++ science).\n"
            << "  --xds-inp then --xinfo then --xparm fill CONTROL before script.\n"
            << "  MOSFLM_IMAGE / CCP4 used when IMAGE/TEMPLATE unset.\n"
            << "  F77 ipmosflm oracle: make fortran-oracle\n";
        return 0;
    }

    // CONTROL input is read up front so the log can start with a copy of it.
    const bool interactive = !script && ::isatty(STDIN_FILENO) != 0;
    std::string control_text;
    if (script) {
        std::ifstream in(*script);
        std::ostringstream text;
        text << in.rdbuf();
        control_text = text.str();
    } else if (!interactive) {
        std::ostringstream text;
        text << std::cin.rdbuf();
        control_text = text.str();
    }

    // --output-dir: a fresh directory for every output, as for ``run``.
    if (output_dir) {
        try {
            output_dir = mosflm::io::create_output_directory(
                output_dir, std::filesystem::current_path(), std::time(nullptr)
            );
        } catch (const std::exception& ex) {
            std::cerr << "ProcessingLab: " << ex.what() << '\n';
            return 2;
        }
        if (log_path && log_path->is_relative()) {
            log_path = *output_dir / *log_path;
        }
    }

    // Log file: everything on stdout/stderr is copied (process.log).
    std::optional<ProcessLog> plog;
    if (log_path) {
        try {
            plog.emplace(
                *log_path,
                interactive ? std::string() : control_input_block(control_text)
            );
        } catch (const std::exception& ex) {
            std::cerr << "ProcessingLab: " << ex.what()
                      << " (continuing without a log file)\n";
        }
    }

    print_banner();
    std::cout << " Run mode:        " << mosflm::io::run_mode_name(run_mode) << "  (parameters from "
              << mosflm::io::run_mode_parameter_source(run_mode) << ")\n";

    std::optional<mosflm::io::ProcessRecord> record;
    try {
        mosflm::control::ControlSession session;
        if (output_dir) {
            session.set_output_directory(*output_dir);
        }
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

        // With --output-dir a script run gets the same records as ``run``:
        // input.json (header, imports, --xparm and CONTROL, per field),
        // process.json with stage lifecycle, process.md.
        if (output_dir && !interactive) {
            try {
                std::string text = control_text;
                if (script) {
                    std::ifstream f(*script);
                    std::ostringstream t;
                    t << f.rdbuf();
                    text = t.str();
                }
                const auto with_imports = [&](mosflm::control::ControlSession& d, bool xparm) {
                    d.set_output_directory(*output_dir);
                    if (xds_inp_path) {
                        mosflm::io::apply_imported_geometry(d, mosflm::io::XdsInpImporter::parse_file(*xds_inp_path));
                    }
                    if (xinfo_path) {
                        mosflm::io::apply_imported_geometry(d, mosflm::io::XinfoImporter::parse_file(*xinfo_path));
                    }
                    if (xparm && xparm_path) {
                        mosflm::io::apply_xparm(d, mosflm::io::read_xparm(*xparm_path));
                    }
                };
                mosflm::control::ControlSession imports_only;
                with_imports(imports_only, true);
                mosflm::control::ControlSession dry;
                with_imports(dry, true);
                dry_parse(dry, text);
                mosflm::io::ExperimentInputs in;
                in.image = img_override ? *img_override
                                        : dry.resolve_image_path(std::nullopt, dry.state().process_first);
                in.xds_inp = xds_inp_path;
                in.xinfo = xinfo_path;
                auto model = mosflm::io::resolve_experiment(in);
                mosflm::io::apply_control_state(model, imports_only.state(), mosflm::io::MetadataSource::Xparm);
                mosflm::io::apply_control_state(model, dry.state(), mosflm::io::MetadataSource::Control);
                std::vector<std::string> cl;
                for (int i = 0; i < argc; ++i) {
                    cl.emplace_back(argv[i]);
                }
                const mosflm::io::RunContext ctx{"script", run_mode, cl, text, xds_inp_path, xinfo_path,
                                                 xparm_path, std::nullopt};
                record.emplace(*output_dir, mosflm::io::input_record(model, dry.state(), ctx),
                               mosflm::io::planned_stages(text));
                observe(session, *record);
            } catch (const std::exception& ex) {
                std::cerr << " *** WARNING **   input.json / process.json not written: " << ex.what() << '\n';
            }
        }

        std::vector<mosflm::control::ControlGoResult> runs;
        if (script) {
            runs = session.run_script(*script, img_override);
        } else {
            // Interactive terminal: say we are waiting and prompt per line.
            if (interactive) {
                std::cout << "\n Reading CONTROL commands from the terminal."
                          << "\n Enter keywords (e.g. FINDSPOTS ..., AUTOINDEX"
                          << " ...), GO to run, EXIT or Ctrl-D to finish.\n\n";
            }
            if (interactive) {
                runs = session.run_stream(std::cin, img_override, " ProcessingLab> ");
            } else {
                // Piped / here-doc input (already read for the log).
                std::istringstream in(control_text);
                runs = session.run_stream(in, img_override, "");
            }
            if (interactive && !session.exit_requested()) {
                std::cout << '\n';  // Ctrl-D leaves the cursor after the prompt.
            }
        }

        std::filesystem::path out_mtz = hklout ? *hklout : session.state().hklout;
        if (output_dir && out_mtz.is_relative()) {
            out_mtz = *output_dir / out_mtz;
        }
        report_runs(session, runs, out_mtz);
        if (record) {
            record->add_output(out_mtz, "mtz");
            record->add_output(*output_dir / "index.log", "index_summary");
            if (log_path) {
                record->add_output(*log_path, "log");
            }
            record->complete();
        }
        if (output_dir) {
            std::cout << "  Output directory: " << output_dir->string() << '\n';
        }

        if (session.needs_no_run_warning()) {
            print_no_run_and_end();
        } else {
            print_end();
        }
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "ProcessingLab: " << ex.what() << '\n';
        if (record) {
            record->fail(ex.what());
        }
        return 1;
    }
}
