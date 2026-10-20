#include "mosflm/control/ControlSession.hpp"

#include "mosflm/crystal/SpaceGroupInfo.hpp"
#include "mosflm/diffraction/DataQualityStats.hpp"
#include "mosflm/diffraction/DetectorInit.hpp"
#include "mosflm/diffraction/TrustedRegion.hpp"
#include "mosflm/geometry/CellFix.hpp"
#include "mosflm/geometry/LabFrame.hpp"
#include "mosflm/io/XparmIO.hpp"
#include "mosflm/image/ImageIO.hpp"
#include "mosflm/refinement/PostRefinement.hpp"
#include "mosflm/util/OpenMP.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mosflm::control {
namespace {

AutoIndexRunner& autoindex_runner_slot() {
    static AutoIndexRunner runner;
    return runner;
}

FindSpotsRunner& findspots_runner_slot() {
    static FindSpotsRunner runner;
    return runner;
}

std::string strip_comment(std::string_view line) {
    // Only ``!`` starts a comment. ``#`` is the MOSFLM TEMPLATE digit
    // placeholder (e.g. ``insulin_1_###.img``) and must be preserved.
    std::string out;
    out.reserve(line.size());
    for (char c : line) {
        if (c == '!') {
            break;
        }
        out.push_back(c);
    }
    while (!out.empty() &&
           std::isspace(static_cast<unsigned char>(out.back()))) {
        out.pop_back();
    }
    return out;
}

std::string upper(std::string s) {
    for (char& c : s) {
        c = static_cast<char>(
            std::toupper(static_cast<unsigned char>(c))
        );
    }
    return s;
}

std::vector<std::string> tokenize(const std::string& line) {
    std::vector<std::string> tok;
    std::istringstream iss(line);
    std::string w;
    while (iss >> w) {
        if (w.size() >= 2 &&
            ((w.front() == '"' && w.back() == '"') ||
             (w.front() == '\'' && w.back() == '\''))) {
            w = w.substr(1, w.size() - 2);
        }
        tok.push_back(w);
    }
    return tok;
}

double require_double(const std::string& s, const char* ctx) {
    try {
        std::size_t idx = 0;
        const double v = std::stod(s, &idx);
        if (idx != s.size()) {
            throw std::invalid_argument("trailing junk");
        }
        return v;
    } catch (const std::exception&) {
        throw std::invalid_argument(
            std::string("CONTROL ") + ctx + ": expected number, got '" + s +
            "'"
        );
    }
}

int require_int(const std::string& s, const char* ctx) {
    try {
        std::size_t idx = 0;
        const int v = std::stoi(s, &idx);
        if (idx != s.size()) {
            throw std::invalid_argument("trailing junk");
        }
        return v;
    } catch (const std::exception&) {
        throw std::invalid_argument(
            std::string("CONTROL ") + ctx + ": expected integer, got '" + s +
            "'"
        );
    }
}

std::string expand_template(const std::string& tmpl, int frame) {
    const auto hash = tmpl.find('#');
    if (hash == std::string::npos) {
        return tmpl;
    }
    std::size_t n = 0;
    while (hash + n < tmpl.size() && tmpl[hash + n] == '#') {
        ++n;
    }
    std::string num = std::to_string(frame);
    if (num.size() < n) {
        num.insert(0, n - num.size(), '0');
    }
    return tmpl.substr(0, hash) + num + tmpl.substr(hash + n);
}

bool is_unimplemented_major(const std::string& key) {
    static const char* k[] = {
        "REFIX", "PICK", "XDISPLAY", "TCPIP", "SOCKET",
        nullptr
    };
    for (const char** p = k; *p != nullptr; ++p) {
        if (key == *p) {
            return true;
        }
    }
    return false;
}

}  // namespace

void ControlSession::register_autoindex_runner(AutoIndexRunner runner) {
    autoindex_runner_slot() = std::move(runner);
}

void ControlSession::register_findspots_runner(FindSpotsRunner runner) {
    findspots_runner_slot() = std::move(runner);
}

void ControlSession::ingest_matrix_numbers(const std::vector<double>& vals) {
    for (double v : vals) {
        if (matrix_needed_ <= 0) {
            break;
        }
        matrix_buf_.push_back(v);
        --matrix_needed_;
    }
    finish_matrix_if_ready();
}

void ControlSession::finish_matrix_if_ready() {
    if (matrix_needed_ != 0 || matrix_buf_.size() != 9) {
        return;
    }
    state_.amat << matrix_buf_[0], matrix_buf_[1], matrix_buf_[2],
        matrix_buf_[3], matrix_buf_[4], matrix_buf_[5], matrix_buf_[6],
        matrix_buf_[7], matrix_buf_[8];
    state_.has_matrix = true;
    matrix_buf_.clear();
}

void ControlSession::parse_autoindex_args(const std::vector<std::string>& args) {
    state_.autoindex_pending = true;
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string u = upper(args[i]);
        if (u == "DPS" || u == "UNITCELL" || u == "UNITCELLDIST" ||
            u == "UNITECELLDIST") {
            continue;
        }
        if ((u == "THRESHOLD" || u == "THRESH") && i + 1 < args.size()) {
            state_.autoindex.min_isigi =
                static_cast<float>(require_double(args[++i], "AUTOINDEX"));
            continue;
        }
        if (u == "HKLDEV" && i + 1 < args.size()) {
            state_.autoindex.hkldev =
                static_cast<float>(require_double(args[++i], "AUTOINDEX"));
            continue;
        }
        if (u == "NVECT" && i + 1 < args.size()) {
            state_.autoindex.nvect = require_int(args[++i], "AUTOINDEX");
            continue;
        }
        if (u == "NVREF" && i + 1 < args.size()) {
            state_.autoindex.nvref = require_int(args[++i], "AUTOINDEX");
            continue;
        }
        if (u == "SOLUTION" && i + 1 < args.size()) {
            // MOSFLM SOLUTION n selects a printed indexing solution; map the
            // common insulin case (12 → cI) via prefer_system. Do not treat n
            // as a Bravais character index (1–44).
            const int sol = require_int(args[++i], "AUTOINDEX");
            if (sol == 12 && state_.autoindex.prefer_system.empty()) {
                state_.autoindex.prefer_system = "cI";
            }
            continue;
        }
        if (u == "IMAGE" && i + 1 < args.size()) {
            state_.autoindex.image_number = require_int(args[++i], "AUTOINDEX");
            continue;
        }
        if (u == "FILE" && i + 1 < args.size()) {
            state_.autoindex.spot_file = args[++i];
            continue;
        }
        if (u == "SYSTEM" && i + 1 < args.size()) {
            state_.autoindex.prefer_system = args[++i];
            continue;
        }
        if (u == "SEPMIN" && i + 1 < args.size()) {
            state_.autoindex.sepmin =
                static_cast<float>(require_double(args[++i], "AUTOINDEX"));
            continue;
        }
        if ((u == "CLUSTER_RADIUS" || u == "CLUSTERRADIUS") &&
            i + 1 < args.size()) {
            state_.autoindex.cluster_radius =
                static_cast<float>(require_double(args[++i], "AUTOINDEX"));
            continue;
        }
        if (u == "MININDEXED" && i + 1 < args.size()) {
            state_.autoindex.min_indexed_frac =
                static_cast<float>(require_double(args[++i], "AUTOINDEX"));
            continue;
        }
        if (u == "DELTA" || u == "DIFFVECTOR" || u == "IDXREF") {
            state_.autoindex.use_diff_vector = true;
            continue;
        }
        if (u == "NODELTA" || u == "NODIFFVECTOR" || u == "DPSONLY") {
            state_.autoindex.use_diff_vector = false;
            continue;
        }
    }
}

void ControlSession::parse_findspots_args(const std::vector<std::string>& args) {
    state_.findspots_pending = true;
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string u = upper(args[i]);
        if (u == "FIND" && i + 1 < args.size()) {
            state_.findspots.image_number = require_int(args[++i], "FINDSPOTS");
            continue;
        }
        if (u == "FILE" && i + 1 < args.size()) {
            state_.findspots.spot_file = args[++i];
            continue;
        }
        if (u == "PHI" && i + 2 < args.size()) {
            state_.findspots.phi_start_deg =
                require_double(args[++i], "FINDSPOTS");
            state_.findspots.phi_end_deg =
                require_double(args[++i], "FINDSPOTS");
            state_.findspots.has_phi = true;
            continue;
        }
        if ((u == "THRESHOLD" || u == "THRESH") && i + 1 < args.size()) {
            state_.findspots.threshold_sigma =
                static_cast<float>(require_double(args[++i], "FINDSPOTS"));
            continue;
        }
        if (u == "MAXSPOTS" && i + 1 < args.size()) {
            state_.findspots.max_spots = require_int(args[++i], "FINDSPOTS");
            continue;
        }
        if (u == "MINDIST" && i + 1 < args.size()) {
            state_.findspots.min_distance = require_int(args[++i], "FINDSPOTS");
            continue;
        }
        if (u == "BORDER" && i + 1 < args.size()) {
            state_.findspots.border = require_int(args[++i], "FINDSPOTS");
            continue;
        }
        if (u == "RANGE" && i + 2 < args.size()) {
            state_.findspots.range_first = require_int(args[++i], "FINDSPOTS");
            state_.findspots.range_last = require_int(args[++i], "FINDSPOTS");
            state_.findspots.multiframe = true;
            continue;
        }
        if (u == "MULTIFRAME" || u == "COLSPOT3D" || u == "COLSPOT") {
            state_.findspots.multiframe = true;
            continue;
        }
        if (u == "MERGE_XY" && i + 1 < args.size()) {
            state_.findspots.merge_xy_mm =
                require_double(args[++i], "FINDSPOTS");
            continue;
        }
        if (u == "MERGE_PHI" && i + 1 < args.size()) {
            state_.findspots.merge_phi_tol_deg =
                require_double(args[++i], "FINDSPOTS");
            continue;
        }
        if (u == "MINFRAMES" && i + 1 < args.size()) {
            state_.findspots.min_frames = require_int(args[++i], "FINDSPOTS");
            continue;
        }
        // Ignore MOSFLM FINDSPOTS noise keywords (LOCAL, AUTORESLN, …).
    }
    // MULTIFRAME without RANGE: default a 50-frame COLSPOT wedge from FIND.
    if (state_.findspots.multiframe && state_.findspots.range_last <= 0) {
        const int start = state_.findspots.image_number > 0
                              ? state_.findspots.image_number
                              : 1;
        state_.findspots.range_first = start;
        state_.findspots.range_last = start + 49;
    }
}

void ControlSession::apply_keyword(
    const std::string& key,
    const std::vector<std::string>& args
) {
    if (key == "DISTANCE") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL DISTANCE needs 1 value");
        }
        state_.distance_mm = require_double(args[0], "DISTANCE");
        state_.has_distance = true;
        return;
    }
    if (key == "WAVELENGTH" || key == "WAVE") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL WAVELENGTH needs 1 value");
        }
        state_.wavelength_A = require_double(args[0], "WAVELENGTH");
        state_.has_wavelength = true;
        return;
    }
    if (key == "PIXEL") {
        if (args.size() < 2) {
            throw std::invalid_argument("CONTROL PIXEL needs 2 values");
        }
        state_.pixel_x_mm = require_double(args[0], "PIXEL");
        state_.pixel_y_mm = require_double(args[1], "PIXEL");
        state_.has_pixel = true;
        return;
    }
    if (key == "BEAM") {
        if (args.size() < 2) {
            throw std::invalid_argument("CONTROL BEAM needs 2 values");
        }
        state_.beam_x_mm = require_double(args[0], "BEAM");
        state_.beam_y_mm = require_double(args[1], "BEAM");
        state_.has_beam = true;
        return;
    }
    if (key == "CELL") {
        if (args.size() < 6) {
            throw std::invalid_argument("CONTROL CELL needs 6 values");
        }
        state_.cell.a = require_double(args[0], "CELL");
        state_.cell.b = require_double(args[1], "CELL");
        state_.cell.c = require_double(args[2], "CELL");
        state_.cell.alpha_deg = require_double(args[3], "CELL");
        state_.cell.beta_deg = require_double(args[4], "CELL");
        state_.cell.gamma_deg = require_double(args[5], "CELL");
        state_.has_cell = true;
        return;
    }
    if (key == "SYMMETRY" || key == "SYMM") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL SYMMETRY needs 1 value");
        }
        try {
            state_.symmetry_number = require_int(args[0], "SYMMETRY");
        } catch (const std::invalid_argument&) {
            const auto* sg =
                crystal::SpaceGroupInfo::find_by_name(args[0].c_str());
            if (sg == nullptr) {
                throw std::invalid_argument(
                    "CONTROL SYMMETRY: unknown space group '" + args[0] + "'"
                );
            }
            state_.symmetry_number = sg->number;
        }
        state_.has_symmetry = true;
        return;
    }
    if (key == "CELLFIX") {
        // CELLFIX [ON|OFF|HARD|SOFT] [softness]
        if (args.empty()) {
            state_.cellfix_enabled = true;
            state_.cellfix_softness = 1.0;
            return;
        }
        const std::string a0 = upper(args[0]);
        if (a0 == "OFF" || a0 == "NO" || a0 == "0") {
            state_.cellfix_enabled = false;
            return;
        }
        state_.cellfix_enabled = true;
        if (a0 == "SOFT") {
            state_.cellfix_softness = 0.5;
        } else if (a0 == "HARD" || a0 == "ON" || a0 == "YES" || a0 == "1") {
            state_.cellfix_softness = 1.0;
        } else {
            state_.cellfix_softness = require_double(args[0], "CELLFIX");
        }
        if (args.size() >= 2) {
            state_.cellfix_softness = require_double(args[1], "CELLFIX");
        }
        return;
    }
    if (key == "MOSAIC" || key == "MOSAICITY") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL MOSAIC needs 1 value");
        }
        state_.mosaic_deg = require_double(args[0], "MOSAIC");
        state_.has_mosaic = true;
        return;
    }
    if (key == "DIVERGENCE" || key == "DIVH" || key == "DIVV") {
        if (key == "DIVERGENCE") {
            if (args.size() < 2) {
                throw std::invalid_argument(
                    "CONTROL DIVERGENCE needs 2 values (DIVH DIVV)"
                );
            }
            state_.divh_deg = require_double(args[0], "DIVERGENCE");
            state_.divv_deg = require_double(args[1], "DIVERGENCE");
        } else if (key == "DIVH") {
            if (args.size() < 1) {
                throw std::invalid_argument("CONTROL DIVH needs 1 value");
            }
            state_.divh_deg = require_double(args[0], "DIVH");
        } else {
            if (args.size() < 1) {
                throw std::invalid_argument("CONTROL DIVV needs 1 value");
            }
            state_.divv_deg = require_double(args[0], "DIVV");
        }
        state_.has_divergence = true;
        return;
    }
    if (key == "SEPARATION") {
        if (args.size() < 2) {
            throw std::invalid_argument("CONTROL SEPARATION needs 2 values");
        }
        state_.separation_x_mm = require_double(args[0], "SEPARATION");
        state_.separation_y_mm = require_double(args[1], "SEPARATION");
        state_.has_separation = true;
        return;
    }
    if (key == "DIRECTORY" || key == "DIR") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL DIRECTORY needs a path");
        }
        state_.directory = args[0];
        return;
    }
    if (key == "TEMPLATE") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL TEMPLATE needs a name");
        }
        state_.image_template = args[0];
        return;
    }
    if (key == "IMAGE") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL IMAGE needs a path");
        }
        state_.image_path = args[0];
        return;
    }
    if (key == "XPARM" || key == "GXPARM") {
        if (args.empty()) {
            throw std::invalid_argument("CONTROL XPARM needs a path");
        }
        const auto g = io::read_xparm(args[0]);
        io::apply_xparm(*this, g);
        std::cout << "  Imported XPARM: " << args[0] << '\n';
        return;
    }
    if (key == "XPARMOUT" || key == "GXPARMOUT") {
        if (args.empty()) {
            throw std::invalid_argument("CONTROL XPARMOUT needs a path");
        }
        const std::string a0 = upper(args[0]);
        if (a0 == "OFF" || a0 == "NO" || a0 == "0") {
            state_.auto_gxparm = false;
            state_.has_xparm_out = false;
            return;
        }
        state_.xparm_out = args[0];
        state_.has_xparm_out = true;
        state_.auto_gxparm = true;
        const auto g = io::xparm_from_control(*this);
        io::write_xparm(state_.xparm_out, g);
        std::cout << "  Wrote XPARM: " << state_.xparm_out << '\n';
        return;
    }
    if (key == "HKLOUT") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL HKLOUT needs a path");
        }
        state_.hklout = args[0];
        return;
    }
    if (key == "MATRIX" || key == "AMAT") {
        matrix_buf_.clear();
        matrix_needed_ = 9;
        std::vector<double> vals;
        vals.reserve(args.size());
        for (const auto& a : args) {
            vals.push_back(require_double(a, key.c_str()));
        }
        ingest_matrix_numbers(vals);
        return;
    }
    if (key == "PROCESS") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL PROCESS needs a frame range");
        }
        state_.process_first = require_int(args[0], "PROCESS");
        state_.process_last = state_.process_first;
        bool saw_last = false;
        for (std::size_t i = 1; i < args.size(); ++i) {
            const std::string u = upper(args[i]);
            if (u == "TO" && i + 1 < args.size()) {
                state_.process_last = require_int(args[i + 1], "PROCESS");
                saw_last = true;
                ++i;
            } else if (u == "START" && i + 1 < args.size()) {
                state_.process_start_phi_deg =
                    require_double(args[i + 1], "PROCESS");
                state_.has_process_phi = true;
                ++i;
            } else if ((u == "ANGLE" || u == "OSCILLATION") &&
                       i + 1 < args.size()) {
                state_.process_angle_deg =
                    require_double(args[i + 1], "PROCESS");
                state_.has_process_phi = true;
                ++i;
            } else if (!saw_last) {
                // MOSFLM shorthand: ``PROCESS 1 100`` ⇒ first=1 last=100.
                try {
                    state_.process_last = require_int(args[i], "PROCESS");
                    saw_last = true;
                } catch (const std::invalid_argument&) {
                    // Non-numeric token; ignore (compat with extra keywords).
                }
            }
        }
        if (state_.process_last < state_.process_first) {
            std::swap(state_.process_last, state_.process_first);
        }
        state_.process_pending = true;
        return;
    }
    if (key == "NMASKS" || key == "PROFILE") {
        // NMASKS nx ny  or  PROFILE NMASKS nx ny
        std::size_t i = 0;
        if (key == "PROFILE" && !args.empty() && upper(args[0]) == "NMASKS") {
            i = 1;
        }
        if (args.size() < i + 2) {
            throw std::invalid_argument("CONTROL NMASKS needs nx ny");
        }
        state_.nmasks_x = require_int(args[i], "NMASKS");
        state_.nmasks_y = require_int(args[i + 1], "NMASKS");
        if (state_.nmasks_x < 1 || state_.nmasks_y < 1) {
            throw std::invalid_argument("CONTROL NMASKS: nx,ny must be >= 1");
        }
        return;
    }
    if (key == "POLARIZATION" || key == "POLARISN") {
        if (args.empty()) {
            throw std::invalid_argument(
                "CONTROL POLARIZATION needs a fraction (0..1)"
            );
        }
        state_.polarization_fraction =
            require_double(args[0], "POLARIZATION");
        return;
    }
    if (key == "AIR" || key == "AIRABS") {
        if (args.empty()) {
            throw std::invalid_argument("CONTROL AIR needs a coefficient");
        }
        state_.air_absorption_mm_inv = require_double(args[0], "AIR");
        return;
    }
    if (key == "PARTIALSCALE" || key == "SCALEPARTIAL") {
        if (args.empty() || upper(args[0]) == "ON" || upper(args[0]) == "YES" ||
            args[0] == "1") {
            state_.scale_by_partiality = true;
        } else {
            state_.scale_by_partiality = false;
        }
        return;
    }
    if (key == "AUTOSIGMA" || key == "AUTOSIGMAS") {
        if (args.empty() || upper(args[0]) == "ON" || upper(args[0]) == "YES" ||
            args[0] == "1") {
            state_.auto_sigma_msd = true;
        } else {
            state_.auto_sigma_msd = false;
        }
        return;
    }
    if (key == "NPROC" || key == "THREADS" || key == "NTHREADS") {
        if (args.empty() || upper(args[0]) == "OFF" || upper(args[0]) == "1") {
            state_.nproc = 1;
        } else {
            state_.nproc = std::max(1, std::min(require_int(args[0], "NPROC"), 32));
        }
        util::apply_nproc(state_.nproc);
        std::cout << "  NPROC: " << state_.nproc << " thread(s)"
                  << (util::openmp_enabled() ? " (OpenMP)" : " (serial build)")
                  << '\n';
        return;
    }
    if (key == "POSTREFINE" || key == "POSTREFL") {
        state_.enable_postrefinement = true;
        for (std::size_t i = 0; i < args.size(); ++i) {
            const std::string u = upper(args[i]);
            if ((u == "PERIOD" || u == "EVERY") && i + 1 < args.size()) {
                state_.postrefine_period =
                    require_int(args[++i], "POSTREFINE");
                continue;
            }
            if (u == "OFF" || u == "NO" || u == "0") {
                state_.enable_postrefinement = false;
                state_.postrefine_period = 0;
            }
            if (u == "ON" || u == "YES") {
                state_.enable_postrefinement = true;
                if (state_.postrefine_period <= 0) {
                    state_.postrefine_period = 25;
                }
            }
        }
        if (!args.empty() && (upper(args[0]) == "OFF" || args[0] == "0")) {
            state_.enable_postrefinement = false;
            state_.postrefine_period = 0;
        }
        return;
    }
    if (key == "NADD") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL NADD needs 1 value");
        }
        state_.nadd = require_int(args[0], "NADD");
        if (state_.nadd < 1) {
            throw std::invalid_argument("CONTROL NADD must be >= 1");
        }
        return;
    }
    if (key == "DMAX") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL DMAX needs 1 value");
        }
        state_.dmax_A = require_double(args[0], "DMAX");
        return;
    }
    if (key == "DMIN") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL DMIN needs 1 value");
        }
        state_.dmin_A = require_double(args[0], "DMIN");
        return;
    }
    if (key == "RESOLUTION") {
        // RESOLUTION dmax dmin  OR  RESOLUTION dmin (high-res only).
        if (args.size() < 1) {
            throw std::invalid_argument(
                "CONTROL RESOLUTION needs 1 or 2 values"
            );
        }
        if (args.size() >= 2) {
            state_.dmax_A = require_double(args[0], "RESOLUTION");
            state_.dmin_A = require_double(args[1], "RESOLUTION");
        } else {
            state_.dmin_A = require_double(args[0], "RESOLUTION");
        }
        return;
    }
    if (key == "YSCALE") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL YSCALE needs 1 value");
        }
        state_.yscale = require_double(args[0], "YSCALE");
        return;
    }
    if (key == "TILT") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL TILT needs 1 value (mrad)");
        }
        state_.tilt_mrad = require_double(args[0], "TILT");
        return;
    }
    if (key == "TWIST") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL TWIST needs 1 value (mrad)");
        }
        state_.twist_mrad = require_double(args[0], "TWIST");
        return;
    }
    if (key == "ROTMAT2") {
        if (args.empty()) {
            state_.use_rotmat2 = true;
            return;
        }
        const std::string u = upper(args[0]);
        if (u == "ON" || u == "YES" || u == "TRUE" || u == "1") {
            state_.use_rotmat2 = true;
        } else if (u == "OFF" || u == "NO" || u == "FALSE" || u == "0") {
            state_.use_rotmat2 = false;
        } else {
            throw std::invalid_argument(
                "CONTROL ROTMAT2 expects ON/OFF (or no argument for ON)"
            );
        }
        return;
    }
    if (key == "MINFRAC" || key == "MINFRACTION" || key == "MIN_FRACTION") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL MINFRAC needs 1 value");
        }
        state_.min_fractioncalc = require_double(args[0], "MINFRAC");
        return;
    }
    if (key == "PHI_OFFSET" || key == "PHIOFFSET") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL PHI_OFFSET needs 1 value");
        }
        state_.phi_offset_deg = require_double(args[0], "PHI_OFFSET");
        return;
    }
    if (key == "PEAKHALF" || key == "PEAK_HALF") {
        if (args.size() < 1) {
            throw std::invalid_argument("CONTROL PEAKHALF needs 1 integer");
        }
        state_.peak_half_px = require_int(args[0], "PEAKHALF");
        if (state_.peak_half_px < 1) {
            throw std::invalid_argument("CONTROL PEAKHALF must be >= 1");
        }
        return;
    }
    if (key == "WIDEN_EWALD" || key == "WIDENEWALD") {
        if (args.empty()) {
            state_.widen_ewald_with_mosaic = true;
            return;
        }
        const std::string u = upper(args[0]);
        if (u == "ON" || u == "YES" || u == "TRUE" || u == "1") {
            state_.widen_ewald_with_mosaic = true;
        } else if (u == "OFF" || u == "NO" || u == "FALSE" || u == "0") {
            state_.widen_ewald_with_mosaic = false;
        } else {
            state_.mosaic_ewald_scale = require_double(args[0], "WIDEN_EWALD");
            state_.widen_ewald_with_mosaic = state_.mosaic_ewald_scale > 0.0;
        }
        return;
    }
    if (key == "AUTOINDEX") {
        parse_autoindex_args(args);
        return;
    }
    if (key == "FINDSPOTS") {
        parse_findspots_args(args);
        return;
    }
    if (key == "STRATEGY") {
        state_.strategy_pending = true;
        for (std::size_t i = 0; i < args.size(); ++i) {
            const std::string u = upper(args[i]);
            if ((u == "ANGLE" || u == "RANGE") && i + 1 < args.size()) {
                state_.strategy_osc_range_deg =
                    require_double(args[++i], "STRATEGY");
            }
        }
        return;
    }
    if (key == "MULTILATTICE" || key == "LATTICE") {
        // MULTILATTICE followed by 9 matrix numbers (extra lattice).
        if (args.size() < 9) {
            throw std::invalid_argument(
                "CONTROL MULTILATTICE needs 9 orientation-matrix values"
            );
        }
        Eigen::Matrix3d m;
        m << require_double(args[0], "MULTILATTICE"),
            require_double(args[1], "MULTILATTICE"),
            require_double(args[2], "MULTILATTICE"),
            require_double(args[3], "MULTILATTICE"),
            require_double(args[4], "MULTILATTICE"),
            require_double(args[5], "MULTILATTICE"),
            require_double(args[6], "MULTILATTICE"),
            require_double(args[7], "MULTILATTICE"),
            require_double(args[8], "MULTILATTICE");
        state_.extra_lattices.push_back(m);
        return;
    }
    if (key == "GO" || key == "RUN") {
        return;
    }
    if (key == "EXIT" || key == "QUIT") {
        state_.exit_requested = true;
        return;
    }
    if (is_unimplemented_major(key)) {
        throw std::runtime_error(
            "CONTROL: keyword '" + key +
            "' is not yet implemented in the C++ core (Phase E)"
        );
    }
    // Harmless / deferred keywords (GENFILE, NEWMAT, TITLE, …): warn once.
    if (warned_unknown_.insert(key).second) {
        std::cerr << " *** WARNING **   CONTROL: ignoring unknown keyword "
                  << key << '\n';
    }
}

void ControlSession::parse_line(std::string_view raw) {
    const std::string line = strip_comment(raw);
    if (line.empty()) {
        return;
    }

    if (matrix_needed_ > 0) {
        const auto tok = tokenize(line);
        std::vector<double> vals;
        vals.reserve(tok.size());
        for (const auto& t : tok) {
            vals.push_back(require_double(t, "MATRIX"));
        }
        ingest_matrix_numbers(vals);
        return;
    }

    const auto tok = tokenize(line);
    if (tok.empty()) {
        return;
    }
    const std::string key = upper(tok[0]);
    std::vector<std::string> args(tok.begin() + 1, tok.end());
    apply_keyword(key, args);
}

void ControlSession::parse_file(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("CONTROL: cannot open " + path.string());
    }
    std::string line;
    while (std::getline(in, line)) {
        parse_line(line);
        if (state_.exit_requested) {
            break;
        }
    }
    if (matrix_needed_ > 0) {
        throw std::runtime_error("CONTROL: incomplete MATRIX (need 9 numbers)");
    }
}

namespace {

bool path_looks_hdf5(const std::filesystem::path& path) {
    std::string e = path.extension().string();
    for (char& c : e) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return e == ".h5" || e == ".nxs" || e == ".hdf5";
}

/// XDS-lite auto σM / σD: infer reflecting range from partial FRACTIONCALC
/// and oscillation width (rocking ≈ osc / frac for partials).
bool auto_estimate_sigma_msd(
    ControlState& state,
    const diffraction::ProcessResult& acc,
    double osc_w_deg
) {
    if (!state.auto_sigma_msd || !(osc_w_deg > 1e-6)) {
        return false;
    }
    std::vector<double> rockings;
    rockings.reserve(256);
    for (const auto& row : acc.mosflm_reflections) {
        if (row.flag != 0 || !(row.sigma > 0.0f)) {
            continue;
        }
        if (static_cast<double>(row.intensity / row.sigma) < 5.0) {
            continue;
        }
        const double f = static_cast<double>(row.fractioncalc);
        if (f > 0.08 && f < 0.95) {
            rockings.push_back(osc_w_deg / f);
        }
    }
    if (rockings.size() < 40) {
        return false;
    }
    const std::size_t mid = rockings.size() / 2U;
    std::nth_element(
        rockings.begin(), rockings.begin() + static_cast<std::ptrdiff_t>(mid),
        rockings.end()
    );
    const double rocking = rockings[mid];
    // Split rocking into mosaic + small divergence (XDS σM / σD stand-in).
    const double div = std::clamp(0.15 * rocking, 0.01, 0.25);
    double mosaic = std::clamp(rocking - div, 0.05, 2.0);
    // If MOSAIC was set in CONTROL, only seed DIV; do not rewrite mosaic.
    if (state.has_mosaic) {
        if (!state.has_divergence) {
            state.divh_deg = div;
            state.divv_deg = div;
            state.has_divergence = true;
            std::cout << std::fixed << std::setprecision(3);
            std::cout << "  AUTOSIGMA: DIV=" << div
                      << " (MOSAIC kept " << state.mosaic_deg << ")\n";
            return true;
        }
        return false;
    }
    const bool changed = std::fabs(mosaic - state.mosaic_deg) > 0.02;
    if (!changed) {
        return false;
    }
    state.mosaic_deg = mosaic;
    state.has_mosaic = true;
    state.divh_deg = div;
    state.divv_deg = div;
    state.has_divergence = true;
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  AUTOSIGMA: MOSAIC=" << mosaic << " DIV=" << div
              << " (from " << rockings.size() << " partials)\n";
    return true;
}

/// CORRECT/INTEGRATE-lite: refine misset from measured XDET/YDET/ROT and
/// write the updated MOSFLM AMAT back into CONTROL state.
bool postrefine_from_observations(
    ControlState& state,
    const diffraction::ProcessResult& one,
    const DetectorGeometry& det
) {
    using refinement::PostRefObservation;
    using refinement::PostRefinementSettings;
    using refinement::PostRefiner;

    std::vector<PostRefObservation> obs;
    obs.reserve(one.mosflm_reflections.size());
    for (const auto& row : one.mosflm_reflections) {
        if (row.flag != 0 || !(row.sigma > 0.0f)) {
            continue;
        }
        const double isigi =
            static_cast<double>(row.intensity / row.sigma);
        if (isigi < 1.5) {
            continue;
        }
        PostRefObservation o;
        o.h = row.h;
        o.k = row.k;
        o.l = row.l;
        o.phi_obs_deg = static_cast<double>(row.rot);
        o.weight = (row.fractioncalc < 0.95f) ? 1.5 : 0.5;
        o.xdet_obs_px = static_cast<double>(row.xdet);
        o.ydet_obs_px = static_cast<double>(row.ydet);
        obs.push_back(o);
    }
    if (obs.size() < 12) {
        return false;
    }

    PostRefinementSettings settings;
    settings.refine_against_detector = true;
    settings.refine_cell_scale = false;  // avoid lattice walk on long wedges
    settings.detector = det;
    settings.beam_x_mm = state.has_beam ? state.beam_x_mm : det.beam_x_mm;
    settings.beam_y_mm = state.has_beam ? state.beam_y_mm : det.beam_y_mm;
    settings.max_iterations = 8;
    // Fine-φ Eiger: predicted ROT ≈ observed; widen φ gate, rely on XDET.
    settings.outlier_phi_deg = 5.0;

    crystal::CrystalMetadata crystal;
    crystal.cell = gemmi::UnitCell(
        state.cell.a, state.cell.b, state.cell.c, state.cell.alpha_deg,
        state.cell.beta_deg, state.cell.gamma_deg
    );
    if (state.has_symmetry) {
        crystal.spacegroup =
            crystal::SpaceGroupInfo::find_by_number(state.symmetry_number);
    }
    const auto orient =
        OrientationMatrix::from_mosflm_orientation_matrix(state.amat);
    const auto result = PostRefiner{}.refine(
        crystal, orient, state.wavelength_A, state.mosaic_deg, obs, settings
    );
    if (!result.ok) {
        return false;
    }
    // Only accept if positional residual improves (or φ when no pos term).
    const bool pos_ok =
        result.rms_pos_initial_px <= 0.0 ||
        result.rms_pos_px <= result.rms_pos_initial_px * 1.02;
    const bool phi_ok =
        result.rms_phi_deg <= result.rms_phi_initial_deg * 1.02;
    const bool accept = pos_ok && phi_ok;
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  POSTREFINE: n=" << result.n_used
              << " rms_phi=" << result.rms_phi_initial_deg << "→"
              << result.rms_phi_deg << " rms_pos="
              << result.rms_pos_initial_px << "→" << result.rms_pos_px
              << " px" << (accept ? " (applied)" : " (rejected)") << '\n';
    if (!accept) {
        return false;
    }
    // Cap misset so a bad batch cannot walk the orientation off lattice.
    const double cap = 0.35;  // degrees
    const auto& p = result.parameters;
    if (std::fabs(p.dphix_deg) > cap || std::fabs(p.dphiy_deg) > cap ||
        std::fabs(p.dphiz_deg) > cap) {
        std::cout << "  POSTREFINE: misset exceeds " << cap
                  << " deg — not applied\n";
        return false;
    }
    state.amat =
        cpp_ub_to_mosflm_amat(result.orientation.ub(), state.wavelength_A);
    state.has_matrix = true;
    // Do not overwrite CONTROL MOSAIC from postref (φ residuals ≠ σM).
    return true;
}

void append_process_result(
    diffraction::ProcessResult& acc,
    diffraction::ProcessResult&& one,
    int batch
) {
    for (auto& row : one.mosflm_reflections) {
        row.batch = batch;
    }
    acc.n_overlapped += one.n_overlapped;
    acc.n_profile_ok += one.n_profile_ok;
    acc.n_omitted += one.n_omitted;
    acc.predictions.insert(
        acc.predictions.end(),
        std::make_move_iterator(one.predictions.begin()),
        std::make_move_iterator(one.predictions.end())
    );
    acc.integrated.insert(
        acc.integrated.end(),
        std::make_move_iterator(one.integrated.begin()),
        std::make_move_iterator(one.integrated.end())
    );
    acc.meas_observations.insert(
        acc.meas_observations.end(),
        std::make_move_iterator(one.meas_observations.begin()),
        std::make_move_iterator(one.meas_observations.end())
    );
    acc.profile_observations.insert(
        acc.profile_observations.end(),
        std::make_move_iterator(one.profile_observations.begin()),
        std::make_move_iterator(one.profile_observations.end())
    );
    acc.mtz_reflections.insert(
        acc.mtz_reflections.end(),
        std::make_move_iterator(one.mtz_reflections.begin()),
        std::make_move_iterator(one.mtz_reflections.end())
    );
    acc.mosflm_reflections.insert(
        acc.mosflm_reflections.end(),
        std::make_move_iterator(one.mosflm_reflections.begin()),
        std::make_move_iterator(one.mosflm_reflections.end())
    );
    // Keep last frame's postrefinement summary if any.
    if (one.postrefinement.ok) {
        acc.postrefinement = std::move(one.postrefinement);
    }
}

/// Per-image counters for XDS INTEGRATE-style PROCESS logging.
struct ProcessImageStats {
    int ier{0};
    double scale{1.0};
    std::size_t nbkg{0};
    std::size_t novl{0};
    std::size_t newald{0};
    std::size_t nstrong{0};
    std::size_t nrej{0};
    double sigmab{0.0};
    double sigmar{0.0};
    double mean_isigi{0.0};
    std::size_t n_integrated{0};
};

ProcessImageStats make_process_image_stats(
    const diffraction::ProcessResult& one,
    double sigmab,
    double sigmar
) {
    ProcessImageStats s;
    s.sigmab = sigmab;
    s.sigmar = sigmar;
    s.newald = one.predictions.size();
    s.novl = one.n_overlapped;
    s.nrej = one.n_omitted;
    s.nstrong = one.n_profile_ok;
    s.n_integrated = one.mtz_reflections.size();

    for (const auto& obs : one.integrated) {
        if (obs.npix_bg > 0) {
            s.nbkg += static_cast<std::size_t>(obs.npix_bg);
        }
    }
    if (s.nbkg == 0) {
        for (const auto& obs : one.meas_observations) {
            if (obs.npix_bg > 0) {
                s.nbkg += static_cast<std::size_t>(obs.npix_bg);
            }
        }
    }

    double sum_isigi = 0.0;
    std::size_t n_isigi = 0;
    for (const auto& row : one.mosflm_reflections) {
        if (row.flag != 0) {
            continue;
        }
        if (!(row.sigma > 0.0f) || !std::isfinite(row.sigma) ||
            !std::isfinite(row.intensity)) {
            continue;
        }
        sum_isigi += static_cast<double>(row.intensity / row.sigma);
        ++n_isigi;
    }
    if (n_isigi > 0) {
        s.mean_isigi = sum_isigi / static_cast<double>(n_isigi);
    }
    if (s.nstrong == 0 && n_isigi > 0) {
        // Fallback when profile learning did not fill n_profile_ok.
        for (const auto& row : one.mosflm_reflections) {
            if (row.flag == 0 && row.sigma > 0.0f &&
                row.intensity / row.sigma >= 3.0f) {
                ++s.nstrong;
            }
        }
    }
    return s;
}

void print_process_integrate_header(
    const ControlState& state,
    int first,
    int last,
    const std::filesystem::path& image_example,
    double osc_range_deg,
    double start_angle_deg,
    const image::ImageMetadata& img_meta
) {
    // Prefer image-header geometry when CONTROL did not override.
    const double wavelength =
        state.has_wavelength ? state.wavelength_A
                             : (img_meta.wavelength_angstrom > 0.0
                                    ? img_meta.wavelength_angstrom
                                    : state.wavelength_A);
    const double distance =
        state.has_distance ? state.distance_mm
                           : (img_meta.detector_distance_mm > 0.0
                                  ? img_meta.detector_distance_mm
                                  : state.distance_mm);
    const double pixel_x =
        state.has_pixel ? state.pixel_x_mm
                        : (img_meta.pixel_size_x_mm > 0.0
                               ? img_meta.pixel_size_x_mm
                               : state.pixel_x_mm);
    const double pixel_y =
        state.has_pixel ? state.pixel_y_mm
                        : (img_meta.pixel_size_y_mm > 0.0
                               ? img_meta.pixel_size_y_mm
                               : state.pixel_y_mm);
    const double beam_x =
        state.has_beam ? state.beam_x_mm
                       : (img_meta.beam_x_mm > 0.0 ? img_meta.beam_x_mm
                                                   : state.beam_x_mm);
    const double beam_y =
        state.has_beam ? state.beam_y_mm
                       : (img_meta.beam_y_mm > 0.0 ? img_meta.beam_y_mm
                                                   : state.beam_y_mm);

    // XDS INTEGRATE-style preamble (honest ProcessingLab field mapping).
    std::cout << "\n ***** INTEGRATE ***** (ProcessingLab PROCESS)\n\n";
    std::cout << " NAME_TEMPLATE_OF_DATA_FRAMES=" << image_example.string()
              << "\n";
    std::cout << " DATA_RANGE= " << std::setw(7) << first << " " << std::setw(7)
              << last << "\n";
    std::cout << " PROFILE_FITTING= TRUE\n";
    std::cout << " SPACE_GROUP_NUMBER= " << std::setw(4)
              << state.symmetry_number << "\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << " UNIT_CELL_CONSTANTS= " << std::setw(8) << state.cell.a << " "
              << std::setw(8) << state.cell.b << " " << std::setw(8)
              << state.cell.c << " " << std::setw(7) << state.cell.alpha_deg
              << " " << std::setw(7) << state.cell.beta_deg << " " << std::setw(7)
              << state.cell.gamma_deg << "\n";
    std::cout << std::setprecision(6);
    std::cout << " OSCILLATION_RANGE= " << std::setw(10) << osc_range_deg
              << " DEGREES\n";
    std::cout << std::setprecision(3);
    std::cout << " STARTING_ANGLE= " << std::setw(8) << start_angle_deg
              << " STARTING_FRAME= " << std::setw(7) << first << "\n";
    std::cout << std::setprecision(5);
    std::cout << " X-RAY_WAVELENGTH= " << std::setw(9) << wavelength
              << " ANGSTROM\n";
    std::cout << std::setprecision(3);
    std::cout << " DETECTOR_DISTANCE= " << std::setw(9) << distance << "\n";
    const double orgx = (pixel_x > 0.0) ? (beam_x / pixel_x) : 0.0;
    const double orgy = (pixel_y > 0.0) ? (beam_y / pixel_y) : 0.0;
    std::cout << " ORGX= " << std::setw(9) << orgx << "  ORGY= " << std::setw(9)
              << orgy << "  (beam / pixel)\n";
    std::cout << std::setprecision(5);
    std::cout << " QX= " << std::setw(8) << pixel_x << " QY= " << std::setw(8)
              << pixel_y << " (mm)\n";
    std::cout << std::setprecision(3);
    std::cout << " BEAM_DIVERGENCE_E.S.D.= ";
    if (state.has_divergence) {
        std::cout << std::setw(8)
                  << (0.5 * std::max(state.divh_deg, state.divv_deg));
        std::cout << "  (from DIVERGENCE " << std::setw(6) << state.divh_deg
                  << " " << std::setw(6) << state.divv_deg << ")\n";
    } else {
        std::cout << std::setw(8) << 0.0 << "  (DIVERGENCE not set)\n";
    }
    std::cout << " REFLECTING_RANGE_E.S.D.= " << std::setw(8) << state.mosaic_deg
              << "  (MOSAICITY)\n";
    std::cout << " RESOLUTION RANGE (ANGSTROM) " << std::setw(9) << state.dmax_A
              << " " << std::setw(9) << state.dmin_A << "\n";
    std::cout << " NMASKS= " << state.nmasks_x << " " << state.nmasks_y << "\n";
    std::cout << " YSCALE= " << std::setw(8) << state.yscale
              << "  TILT= " << std::setw(8) << state.tilt_mrad
              << "  TWIST= " << std::setw(8) << state.twist_mrad
              << "  WIDEN_EWALD= " << (state.widen_ewald_with_mosaic ? "ON" : "OFF")
              << "\n";
    std::cout << std::flush;
}

void print_process_batch_banner(
    int batch_first,
    int batch_last,
    bool print_symbols
) {
    std::cout << "\n ******************************************************************************\n";
    std::cout << "                     PROCESSING OF IMAGES " << std::setw(7)
              << batch_first << " ... " << std::setw(7) << batch_last << "\n";
    std::cout << " ******************************************************************************\n\n";
    std::cout << " USING   1 PROCESSORS\n\n";
    if (print_symbols) {
        // Column meanings follow XDS INTEGRATE; MEANIS is ProcessingLab-only.
        std::cout << " *** DEFINITION OF SYMBOLS ***\n";
        std::cout << " IER     = ERROR CODE AFTER ACCESSING DATA IMAGE\n";
        std::cout << "              0: NO ERROR\n";
        std::cout << "             -1: CANNOT OPEN OR READ IMAGE FILE\n";
        std::cout << "             -3: WRONG DATA FORMAT / UNSUPPORTED IMAGE\n";
        std::cout << " SCALE   = SCALING FACTOR FOR THIS DATA IMAGE\n";
        std::cout << "            (1.0 until SCALE / postrefinement)\n";
        std::cout << " NBKG    = SUM OF BACKGROUND PIXELS USED PER REFLECTION\n";
        std::cout << " NOVL    = NUMBER OF OVERLAPPED REFLECTIONS ON DATA IMAGE\n";
        std::cout << "            (XDS: overloaded; here: spatial overlap)\n";
        std::cout << " NEWALD  = NUMBER OF PREDICTED REFLECTIONS (on detector)\n";
        std::cout << " NSTRONG = NUMBER OF PROFILE-FITTED / STRONG REFLECTIONS\n";
        std::cout << " NREJ    = NUMBER OF OMITTED REFLECTIONS\n";
        std::cout << " SIGMAB  = BEAM_DIVERGENCE_E.S.D.=SIGMAB (input; deg)\n";
        std::cout << " SIGMAR  = REFLECTING_RANGE_E.S.D.=SIGMAR (MOSAICITY)\n";
        std::cout << " MEANIS  = MEAN I/SIGMA OF ACCEPTED REFLECTIONS ON IMAGE\n";
        std::cout << "            (ProcessingLab extension; not an XDS column)\n\n";
    }
    std::cout << " IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ  SIGMAB  SIGMAR  MEANIS\n";
    std::cout << std::flush;
}

void print_process_image_row(int frame, const ProcessImageStats& s) {
    std::cout << std::setw(6) << frame << " " << std::setw(3) << s.ier << " "
              << std::fixed << std::setprecision(3) << std::setw(6) << s.scale
              << " " << std::setw(7) << s.nbkg << " " << std::setw(4) << s.novl
              << " " << std::setw(6) << s.newald << " " << std::setw(7)
              << s.nstrong << " " << std::setw(5) << s.nrej << " "
              << std::setprecision(4) << std::setw(7) << s.sigmab << " "
              << std::setw(7) << s.sigmar << " " << std::setprecision(1)
              << std::setw(7) << s.mean_isigi << "\n"
              << std::flush;
}

void print_process_batch_footer(
    std::size_t n_integrated_batch,
    const ControlState& state
) {
    std::cout << "\n  " << n_integrated_batch
              << " REFLECTION INTENSITIES INTEGRATED\n";
    std::cout << " SPACE GROUP NUMBER " << std::setw(6) << state.symmetry_number
              << "\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << " UNIT CELL PARAMETERS " << std::setw(8) << state.cell.a << " "
              << std::setw(8) << state.cell.b << " " << std::setw(8)
              << state.cell.c << " " << std::setw(7) << state.cell.alpha_deg
              << " " << std::setw(7) << state.cell.beta_deg << " " << std::setw(7)
              << state.cell.gamma_deg << "\n";
    std::cout << " CRYSTAL MOSAICITY (DEGREES) " << std::setw(8)
              << state.mosaic_deg << "\n";
    std::cout << " CRYSTAL TO DETECTOR DISTANCE (mm) " << std::setw(10)
              << state.distance_mm << "\n";
    std::cout << std::flush;
}

void print_process_run_summary(
    const diffraction::ProcessResult& acc,
    int first,
    int last,
    std::size_t n_frames,
    const crystal::CrystalMetadata& crystal,
    double dmax_A,
    double dmin_A
) {
    std::cout << "\n ******************************************************************************\n";
    std::cout << "  PROCESS SUMMARY  FRAMES " << first << " ... " << last << "  ("
              << n_frames << " image(s))\n";
    std::cout << " ******************************************************************************\n";
    std::cout << "  " << std::setw(8) << acc.mtz_reflections.size()
              << " REFLECTIONS INTEGRATED\n";
    std::cout << "  " << std::setw(8) << acc.predictions.size()
              << " PREDICTIONS (all frames, concatenated)\n";
    std::cout << "  " << std::setw(8) << acc.n_profile_ok
              << " PROFILE-FITTED (strong)\n";
    std::cout << "  " << std::setw(8) << acc.n_overlapped << " OVERLAPPED\n";
    std::cout << "  " << std::setw(8) << acc.n_omitted << " OMITTED\n";
    if (!acc.mosflm_reflections.empty()) {
        double sum = 0.0;
        std::size_t n = 0;
        for (const auto& row : acc.mosflm_reflections) {
            if (row.flag == 0 && row.sigma > 0.0f &&
                std::isfinite(row.intensity) && std::isfinite(row.sigma)) {
                sum += static_cast<double>(row.intensity / row.sigma);
                ++n;
            }
        }
        if (n > 0) {
            std::cout << std::fixed << std::setprecision(2);
            std::cout << "  MEAN I/SIGMA (accepted) = "
                      << (sum / static_cast<double>(n)) << "  (N=" << n
                      << ")\n";
        }
    }
    std::cout << "\n" << std::flush;

    if (!acc.mosflm_reflections.empty() && crystal.spacegroup != nullptr) {
        const auto quality = diffraction::compute_data_quality(
            crystal, acc.mosflm_reflections, dmax_A, dmin_A
        );
        diffraction::write_data_quality_table(std::cout, quality);
    }
}

}  // namespace

std::filesystem::path ControlSession::resolve_image_path(
    const std::optional<std::filesystem::path>& image_override,
    int frame
) const {
    const int use_frame = (frame > 0) ? frame : state_.process_first;
    // TEMPLATE must win over --image / IMAGE for multi-frame CBF/SMV: those
    // formats ignore the frame argument to read_diffraction_image, so reusing
    // a single --image path clones frame-1 pixels, φ, XDET and I onto every
    // batch (I/σ collapses to noise about the same empty boxes).
    if (!state_.image_template.empty()) {
        const std::string name =
            expand_template(state_.image_template, use_frame);
        return state_.directory / name;
    }
    if (image_override && !image_override->empty()) {
        return *image_override;
    }
    if (!state_.image_path.empty()) {
        if (state_.image_path.is_absolute()) {
            return state_.image_path;
        }
        return state_.directory / state_.image_path;
    }
    (void)use_frame;
    throw std::runtime_error(
        "CONTROL: no IMAGE/TEMPLATE and no image override for PROCESS"
    );
}

diffraction::ProcessRequest ControlSession::make_process_request(
    const image::DiffractionImage& image,
    int frame
) const {
    if (!state_.has_matrix) {
        throw std::runtime_error(
            "CONTROL: MATRIX/AMAT required before PROCESS GO"
        );
    }
    if (!state_.has_cell) {
        throw std::runtime_error("CONTROL: CELL required before PROCESS GO");
    }

    const gemmi::SpaceGroup* sg = nullptr;
    if (state_.has_symmetry) {
        sg = crystal::SpaceGroupInfo::find_by_number(state_.symmetry_number);
        if (sg == nullptr) {
            throw std::runtime_error(
                "CONTROL: unknown SYMMETRY number " +
                std::to_string(state_.symmetry_number)
            );
        }
    } else {
        sg = crystal::SpaceGroupInfo::find_by_name("P 1");
    }

    // CELLFIX is applied in go_process (mutates state_.cell once). Here we
    // just use the current cell (already fixed if go_process ran CELLFIX).
    crystal::CrystalMetadata crystal;
    crystal.cell = gemmi::UnitCell(
        state_.cell.a,
        state_.cell.b,
        state_.cell.c,
        state_.cell.alpha_deg,
        state_.cell.beta_deg,
        state_.cell.gamma_deg
    );
    if (!crystal.cell.is_crystal()) {
        throw std::runtime_error("CONTROL: invalid CELL");
    }
    crystal.spacegroup = sg;

    DetectorGeometry det = image.metadata.to_detector_geometry();
    det.beam_x_mm = 0.0;
    det.beam_y_mm = 0.0;
    if (state_.has_distance) {
        det.distance_mm = state_.distance_mm;
    } else if (det.distance_mm <= 0.0) {
        det.distance_mm = state_.distance_mm;
    }
    if (state_.has_wavelength) {
        det.wavelength_A = state_.wavelength_A;
    } else if (det.wavelength_A <= 0.0) {
        det.wavelength_A = state_.wavelength_A;
    }
    if (state_.has_pixel) {
        det.pixel_size_x_mm = state_.pixel_x_mm;
        det.pixel_size_y_mm = state_.pixel_y_mm;
    }
    // Detector fiducials: neutral defaults unless CONTROL YSCALE/TILT/TWIST.
    det.yscale = state_.yscale;
    det.tilt_mrad = state_.tilt_mrad;
    det.twist_mrad = state_.twist_mrad;
    det.use_rotmat2 = state_.use_rotmat2;

    diffraction::ScanRange scan;
    if (state_.has_process_phi) {
        const int offset = frame - state_.process_first;
        scan.phi_start_deg =
            state_.process_start_phi_deg +
            static_cast<double>(offset) * state_.process_angle_deg;
        scan.phi_end_deg = scan.phi_start_deg + state_.process_angle_deg;
    } else if (image.metadata.oscillation_range_deg > 0.0) {
        scan.phi_start_deg = image.metadata.oscillation_start_deg;
        scan.phi_end_deg =
            scan.phi_start_deg + image.metadata.oscillation_range_deg;
    } else {
        const int offset = frame - state_.process_first;
        scan.phi_start_deg =
            state_.process_start_phi_deg +
            static_cast<double>(offset) * state_.process_angle_deg;
        scan.phi_end_deg = scan.phi_start_deg + state_.process_angle_deg;
    }
    if (std::fabs(state_.phi_offset_deg) > 1e-12) {
        scan.phi_start_deg += state_.phi_offset_deg;
        scan.phi_end_deg += state_.phi_offset_deg;
    }

    diffraction::PredictionSettings prediction;
    prediction.mosaic_spread_deg = state_.mosaic_deg;
    prediction.divh_deg = state_.has_divergence ? state_.divh_deg : 0.0;
    prediction.divv_deg = state_.has_divergence ? state_.divv_deg : 0.0;
    prediction.dmin_A = state_.dmin_A;
    prediction.dmax_A = state_.dmax_A;
    prediction.phi_step_deg = state_.phi_step_deg;
    prediction.apply_systematic_absences = true;
    // REEK PROCESS: do not widen Ewald with mosaic (avoids multi-frame
    // over-prediction / extreme multiplicity vs DIALS).
    prediction.exclude_ice_rings = false;
    prediction.exclude_cusp = true;
    prediction.cusp_fraction = 0.035;
    // NWMAX-like cap (degrees). Multiplicity is controlled by not widening
    // the Ewald gate with mosaic (REEK PROCESS policy).
    prediction.max_phi_width_deg = 5.0;
    prediction.apply_reek_detector_margin = true;
    prediction.detector_margin_px = 24.0;
    prediction.widen_ewald_with_mosaic = state_.widen_ewald_with_mosaic;
    prediction.mosaic_ewald_scale = state_.mosaic_ewald_scale;
    // Keep Ewald angular half-width well below the image osc so mid-frame
    // spots are classified ``full`` (FRACTIONCALC=1) instead of spanning the
    // whole image as 0.5-partials. Default 0.0025 Å⁻¹ ≈ 0.14° at λ≈0.98 —
    // comparable to a 0.15° frame and made almost everything partial.
    {
        constexpr double kDeg2Rad = 3.14159265358979323846 / 180.0;
        const double osc_w = std::max(
            1e-6, scan.phi_end_deg - scan.phi_start_deg
        );
        // Tuned for thaumatin/DIALS-like obs/image (~500–700) with correct beam.
        const double half_deg = 0.18 * osc_w;  // ~0.027° for 0.15° frames
        const double lam =
            det.wavelength_A > 0.0 ? det.wavelength_A : state_.wavelength_A;
        if (lam > 0.0) {
            prediction.ewald_tolerance_Ainv =
                (half_deg * kDeg2Rad) / lam;
        }
    }

    diffraction::IntegrationSettings integration;
    integration.box_half_px = 3;
    integration.bg_frame_px = 2;
    integration.apply_lp = true;
    integration.polarization_fraction = state_.polarization_fraction;
    integration.mosflm_detector_axes = true;
    integration.image_beam_x_mm =
        state_.has_beam ? state_.beam_x_mm : image.metadata.beam_x_mm;
    integration.image_beam_y_mm =
        state_.has_beam ? state_.beam_y_mm : image.metadata.beam_y_mm;

    diffraction::MeasurementBoxSettings measurement;
    measurement.separation_x_mm = state_.separation_x_mm;
    measurement.separation_y_mm = state_.separation_y_mm;
    // Fine-pixel (Eiger ≤0.08 mm): thicker bg rim; else modest rim.
    const double px_mm = state_.has_pixel
                             ? state_.pixel_x_mm
                             : image.metadata.pixel_size_x_mm;
    // Clamp PEAKHALF so the peak box stays within half SEPARATION (avoids
    // neighbour bleed on fine-pixel detectors).
    {
        const double px = std::max(px_mm, 1e-6);
        // Peak box must fit inside SEPARATION so neighbours do not share
        // peak pixels (half-separation in px is the hard upper bound).
        const double half_sep_px = 0.5 * std::min(
            state_.separation_x_mm, state_.separation_y_mm
        ) / px;
        const int max_half = std::max(
            1, static_cast<int>(std::floor(half_sep_px + 1e-9))
        );
        int ph = state_.peak_half_px >= 1 ? state_.peak_half_px : -1;
        if (ph < 0) {
            ph = std::min(5, max_half);
        }
        measurement.peak_half_px = std::min(ph, max_half);
    }
    // Eiger: modest rim (4); oversized rims dilute SNR without helping CC.
    measurement.bg_rim_px = (px_mm > 0.0 && px_mm <= 0.08) ? 4 : 3;
    measurement.bg_reject_nsigma = 3.0;
    // INIT-lite BLANK/GAIN: tighten bg reject when gain is high (noisy PAD).
    {
        const auto init = diffraction::estimate_blank_gain(image, 32);
        if (init.ok && init.gain > 1.5) {
            measurement.bg_reject_nsigma =
                std::max(2.5, 3.0 / std::sqrt(init.gain));
        }
    }
    measurement.mark_overlaps = true;
    measurement.integrate_overlaps = true;
    // Local COM refine; cap shift to avoid noise walks on weak Eiger spots.
    measurement.refine_centroid = true;
    measurement.centroid_half_px = 2;
    measurement.centroid_max_shift_px = 1.5;
    // INIT/DEFPIX-lite trusted region (circular + DMIN/DMAX).
    {
        diffraction::TrustedRegionSettings tr;
        tr.trusted_radius_frac = 0.98;
        tr.dmin_A = state_.dmin_A;
        tr.dmax_A = state_.dmax_A > 0.0 ? state_.dmax_A : 0.0;
        if (state_.has_beam) {
            tr.beam_x_mm = state_.beam_x_mm;
            tr.beam_y_mm = state_.beam_y_mm;
        }
        measurement.trusted_mask =
            diffraction::build_trusted_mask(image, det, tr);
        measurement.trusted_width = image.metadata.width;
    }

    diffraction::StandardProfileSettings profile;
    profile.box = measurement;
    // Eiger: learn from I/σ≥3; prefer profile only when pkratio is sane.
    profile.learn_min_isigi = 3.0;
    profile.min_learners = 8;
    profile.fallback_sigma_px = 1.5;
    profile.prefer_profile = true;
    profile.max_pkratio_for_profile = 2.0;
    profile.n_bins_x = state_.nmasks_x;
    profile.n_bins_y = state_.nmasks_y;

    refinement::PostRefinementSettings postref;
    postref.refine_against_detector = true;
    postref.detector = det;
    postref.beam_x_mm =
        state_.has_beam ? state_.beam_x_mm : image.metadata.beam_x_mm;
    postref.beam_y_mm =
        state_.has_beam ? state_.beam_y_mm : image.metadata.beam_y_mm;
    // Prefer ControlSession periodic (measured XDET) over pipeline pre-pass.
    const bool pipeline_postref =
        state_.enable_postrefinement && state_.postrefine_period <= 0;

    return diffraction::ProcessRequest{
        crystal,
        OrientationMatrix::from_mosflm_orientation_matrix(state_.amat),
        det,
        scan,
        prediction,
        integration,
        measurement,
        profile,
        true,
        true,
        pipeline_postref,
        postref,
        state_.min_fractioncalc,
        state_.scale_by_partiality,
        state_.air_absorption_mm_inv,
        image
    };
}

ControlGoResult ControlSession::go(
    const std::optional<std::filesystem::path>& image_override
) {
    ControlGoResult out;
    if (matrix_needed_ > 0) {
        throw std::runtime_error("CONTROL: incomplete MATRIX before GO");
    }

    if (state_.findspots_pending) {
        if (!findspots_runner_slot()) {
            throw std::runtime_error(
                "CONTROL FINDSPOTS: spot finder not linked into this build"
            );
        }
        findspots_runner_slot()(state_, image_override);
        state_.findspots_pending = false;
        state_.done_run = true;
        out.ran_findspots = true;
        out.image_used = resolve_image_path(image_override);
        out.message = "FINDSPOTS wrote " +
                      state_.findspots.spot_file.string() + " (" +
                      std::to_string(state_.findspots_n_spots) + " spots)";
        return out;
    }

    if (state_.autoindex_pending) {
        if (!autoindex_runner_slot()) {
            throw std::runtime_error(
                "CONTROL AUTOINDEX: DPS indexer not linked into this build"
            );
        }
        // Prefer FILE from prior FINDSPOTS when AUTOINDEX omitted it.
        if (state_.autoindex.spot_file.empty() &&
            !state_.findspots.spot_file.empty()) {
            state_.autoindex.spot_file = state_.findspots.spot_file;
        }
        autoindex_runner_slot()(state_);
        state_.autoindex_pending = false;
        state_.done_run = true;
        out.ran_autoindex = true;
        out.message = "AUTOINDEX complete";
        // IDXREF analogue: write XPARM.XDS after successful index.
        if (state_.has_matrix && state_.has_cell) {
            std::filesystem::path out_xparm =
                state_.directory / "XPARM.XDS";
            try {
                const auto g = io::xparm_from_control(*this);
                io::write_xparm(out_xparm, g);
                std::cout << "  Wrote XPARM: " << out_xparm << '\n';
            } catch (const std::exception& ex) {
                std::cerr << " *** WARNING **   XPARM write failed: "
                          << ex.what() << '\n';
            }
        }
        return out;
    }

    if (state_.strategy_pending) {
        // Minimal STRATEGY: estimate completeness for cubic/random as
        // osc_range / 90 clipped to 1 (placeholder until full port).
        const double frac = state_.strategy_osc_range_deg / 90.0;
        state_.strategy_completeness_est = std::min(1.0, std::max(0.0, frac));
        state_.strategy_pending = false;
        state_.strategy_done = true;
        state_.done_run = true;
        out.ran_strategy = true;
        out.message = "STRATEGY completeness_est=" +
                      std::to_string(state_.strategy_completeness_est);
        return out;
    }

    if (!state_.process_pending) {
        out.message = "GO with no pending PROCESS/FINDSPOTS/AUTOINDEX/STRATEGY";
        return out;
    }

    const int first = state_.process_first;
    int last = state_.process_last;
    if (first < 1) {
        throw std::runtime_error("CONTROL PROCESS: first image must be >= 1");
    }

    const auto path0 = resolve_image_path(image_override, first);
    out.image_used = path0;

    // Multi-frame SMV/CBF needs a TEMPLATE (or one IMAGE per GO). HDF5 masters
    // reuse the same path with a 1-based frame index.
    const bool hdf5 = path_looks_hdf5(path0) ||
                      (image_override && path_looks_hdf5(*image_override));
    if (hdf5) {
        const int nimg = image::num_diffraction_images(path0);
        if (nimg > 0 && last > nimg) {
            std::cerr << " *** WARNING **   CONTROL PROCESS: last image "
                      << last << " clamped to nimages=" << nimg << '\n';
            last = nimg;
        }
    }
    if (last < first) {
        throw std::runtime_error(
            "CONTROL PROCESS: last image < first after nimages clamp"
        );
    }
    out.process_first = first;
    out.process_last = last;

    // MOSFLM CELLFIX: project cell onto crystal-system metric before integrate.
    if (state_.cellfix_enabled && state_.has_symmetry) {
        const auto* sg =
            crystal::SpaceGroupInfo::find_by_number(state_.symmetry_number);
        if (sg != nullptr) {
            const auto chk =
                check_lattice_vs_spacegroup(state_.cell, sg, 2.0);
            std::cout << "  CELLFIX check: " << chk.message << '\n';
            const UnitCell fixed = cellfix_for_spacegroup(
                state_.cell, sg, state_.cellfix_softness
            );
            if (std::fabs(fixed.a - state_.cell.a) > 1e-4 ||
                std::fabs(fixed.b - state_.cell.b) > 1e-4 ||
                std::fabs(fixed.c - state_.cell.c) > 1e-4 ||
                std::fabs(fixed.alpha_deg - state_.cell.alpha_deg) > 1e-3 ||
                std::fabs(fixed.beta_deg - state_.cell.beta_deg) > 1e-3 ||
                std::fabs(fixed.gamma_deg - state_.cell.gamma_deg) > 1e-3) {
                std::cout << "  CELLFIX applied: " << state_.cell.describe()
                          << " → " << fixed.describe() << '\n';
                state_.cell = fixed;
            }
        }
    }

    const bool multi = last > first;
    if (multi && !hdf5 && state_.image_template.empty()) {
        throw std::runtime_error(
            "CONTROL PROCESS " + std::to_string(first) + " TO " +
            std::to_string(last) +
            ": multi-frame SMV/CBF requires TEMPLATE (HDF5 masters OK)"
        );
    }

    diffraction::ProcessPipeline pipe;
    diffraction::ProcessResult acc;
    crystal::CrystalMetadata crystal;
    std::size_t n_frames = 0;

    // XDS INTEGRATE-style batch size for log sections.
    constexpr int k_log_batch = 25;
    double sigmab =
        state_.has_divergence
            ? 0.5 * std::max(state_.divh_deg, state_.divv_deg)
            : 0.0;
    double sigmar = state_.mosaic_deg;

    // Peek first frame for osc metadata (reuse in the loop).
    const auto image0 = image::read_diffraction_image(path0, first);
    if (!state_.has_wavelength &&
        image0.metadata.wavelength_angstrom > 0.0) {
        state_.wavelength_A = image0.metadata.wavelength_angstrom;
        state_.has_wavelength = true;
    }
    double osc_log = state_.process_angle_deg;
    double start_log = state_.process_start_phi_deg;
    if (!state_.has_process_phi &&
        image0.metadata.oscillation_range_deg > 0.0) {
        osc_log = image0.metadata.oscillation_range_deg;
        start_log = image0.metadata.oscillation_start_deg;
    }
    print_process_integrate_header(
        state_, first, last, path0, osc_log, start_log, image0.metadata
    );

    int batch_first = first;
    std::size_t batch_integrated = 0;
    bool batch_header_printed = false;
    // Multi-frame: periodic postref on by default (period>0); POSTREFINE OFF
    // clears period. Single-frame: only when POSTREFINE explicitly set.
    const bool want_periodic_postref =
        (multi && state_.postrefine_period > 0) ||
        state_.enable_postrefinement;
    bool autosigma_done = false;

    util::apply_nproc(state_.nproc);
    const int nproc = std::max(1, state_.nproc);
    // Chunk size matches POSTREFINE period so AMAT/MOSAIC stay frozen inside
    // an OpenMP region; barriers stay serial.
    const int chunk_size =
        want_periodic_postref
            ? std::max(1, state_.postrefine_period)
            : (last - first + 1);
    if (multi && nproc > 1) {
        std::cout << "  NPROC: PROCESS OpenMP " << nproc
                  << " thread(s), chunk=" << chunk_size << '\n';
    }

    struct FrameOutcome {
        diffraction::ProcessResult result{};
        crystal::CrystalMetadata crystal{};
        ProcessImageStats stats{};
        std::string error;
        bool ok{false};
    };

    for (int chunk0 = first; chunk0 <= last; chunk0 += chunk_size) {
        const int chunk1 = std::min(chunk0 + chunk_size - 1, last);
        const int n_chunk = chunk1 - chunk0 + 1;
        std::vector<FrameOutcome> outcomes(
            static_cast<std::size_t>(n_chunk)
        );

        // Snapshot geometry for this chunk (state_ must not mutate in parallel).
        const double chunk_sigmab = sigmab;
        const double chunk_sigmar = sigmar;
        if (multi && n_chunk > 1) {
            std::cout << "  PROCESS: frames " << chunk0 << "-" << chunk1
                      << " (nproc=" << nproc << ") …" << std::flush;
        }

#if defined(PROCESSINGLAB_HAS_OPENMP) && PROCESSINGLAB_HAS_OPENMP
#pragma omp parallel for num_threads(nproc) schedule(dynamic, 1) \
    if (nproc > 1 && n_chunk > 1)
#endif
        for (int frame = chunk0; frame <= chunk1; ++frame) {
            const std::size_t i =
                static_cast<std::size_t>(frame - chunk0);
            FrameOutcome& outc = outcomes[i];
            try {
                const auto img_path =
                    resolve_image_path(image_override, frame);
                const auto image =
                    (frame == first)
                        ? image0
                        : image::read_diffraction_image(img_path, frame);
                auto req = make_process_request(image, frame);
                outc.crystal = req.crystal;
                diffraction::ProcessPipeline local_pipe;
                auto one = local_pipe.run(req);
                outc.stats =
                    make_process_image_stats(one, chunk_sigmab, chunk_sigmar);
                outc.result = std::move(one);
                outc.ok = true;
            } catch (const std::exception& ex) {
                outc.ok = false;
                outc.error = ex.what();
                outc.stats.ier = -1;
                outc.stats.sigmab = chunk_sigmab;
                outc.stats.sigmar = chunk_sigmar;
            }
        }
        if (multi && n_chunk > 1) {
            std::cout << " done\n";
        }

        // Serial: logging, append, postrefine, AUTOSIGMA (frame order).
        for (int frame = chunk0; frame <= chunk1; ++frame) {
            if (!batch_header_printed ||
                (frame - batch_first) >= k_log_batch) {
                if (batch_header_printed) {
                    print_process_batch_footer(batch_integrated, state_);
                    if (!autosigma_done && state_.auto_sigma_msd) {
                        autosigma_done = auto_estimate_sigma_msd(
                            state_, acc, osc_log
                        );
                        sigmar = state_.mosaic_deg;
                        sigmab = state_.has_divergence
                                     ? 0.5 * std::max(
                                           state_.divh_deg, state_.divv_deg
                                       )
                                     : 0.0;
                    }
                    batch_integrated = 0;
                    batch_first = frame;
                }
                const int batch_last =
                    std::min(batch_first + k_log_batch - 1, last);
                print_process_batch_banner(
                    batch_first, batch_last,
                    /*print_symbols=*/!batch_header_printed
                );
                batch_header_printed = true;
            }

            const std::size_t i =
                static_cast<std::size_t>(frame - chunk0);
            FrameOutcome& outc = outcomes[i];
            if (!outc.ok) {
                print_process_image_row(frame, outc.stats);
                std::cout << " *** IER=-1: " << outc.error << "\n";
                throw std::runtime_error(outc.error);
            }
            batch_integrated += outc.stats.n_integrated;
            if (n_frames == 0) {
                crystal = outc.crystal;
            }
            if (want_periodic_postref) {
                const int period = std::max(1, state_.postrefine_period);
                const bool due =
                    state_.enable_postrefinement
                        ? ((frame - first) % period == 0)
                        : ((frame - first) % period == period - 1 ||
                           frame == last);
                if (due) {
                    const auto img_path =
                        resolve_image_path(image_override, frame);
                    const auto image =
                        (frame == first)
                            ? image0
                            : image::read_diffraction_image(img_path, frame);
                    auto req = make_process_request(image, frame);
                    postrefine_from_observations(
                        state_, outc.result, req.detector
                    );
                }
            }
            print_process_image_row(frame, outc.stats);
            append_process_result(acc, std::move(outc.result), frame);
            ++n_frames;
        }
    }

    if (batch_header_printed) {
        print_process_batch_footer(batch_integrated, state_);
    }
    print_process_run_summary(
        acc, first, last, n_frames, crystal, state_.dmax_A, state_.dmin_A
    );

    out.crystal = crystal;
    out.process = std::move(acc);

    std::string msg = "PROCESS complete frames " + std::to_string(first) +
                      "-" + std::to_string(last) + " (" +
                      std::to_string(n_frames) + " image(s))" +
                      "; integrated=" +
                      std::to_string(out.process.mtz_reflections.size());
    if (state_.nadd > 1 && !state_.image_template.empty() &&
        state_.enable_postrefinement) {
        msg += "; NADD=" + std::to_string(state_.nadd) +
               " reserved for joint refine";
    }
    out.message = std::move(msg);

    // CORRECT-lite: auto-write GXPARM.XDS after PROCESS (XDS GXPARM analogue).
    if (state_.auto_gxparm) {
        std::filesystem::path out_xparm = state_.xparm_out;
        if (out_xparm.empty()) {
            out_xparm = "GXPARM.XDS";
        }
        if (!out_xparm.is_absolute()) {
            out_xparm = state_.directory / out_xparm;
        }
        try {
            const auto g = io::xparm_from_control(*this);
            io::write_xparm(out_xparm, g);
            state_.xparm_out = out_xparm;
            state_.has_xparm_out = true;
            std::cout << "  Wrote GXPARM: " << out_xparm << '\n';
        } catch (const std::exception& ex) {
            std::cerr << " *** WARNING **   GXPARM write failed: " << ex.what()
                      << '\n';
        }
    }

    out.ran_process = true;
    state_.process_pending = false;
    state_.done_run = true;
    return out;
}

std::vector<ControlGoResult> ControlSession::run_stream(
    std::istream& in,
    const std::optional<std::filesystem::path>& image_override
) {
    std::vector<ControlGoResult> results;
    std::string line;
    while (std::getline(in, line)) {
        const std::string stripped = strip_comment(line);
        if (stripped.empty()) {
            continue;
        }
        const auto tok = tokenize(stripped);
        if (tok.empty()) {
            continue;
        }
        const std::string key = upper(tok[0]);
        if (key == "GO" || key == "RUN") {
            if (matrix_needed_ > 0) {
                throw std::runtime_error("CONTROL: incomplete MATRIX before GO");
            }
            results.push_back(go(image_override));
            continue;
        }
        parse_line(line);
        if (state_.exit_requested) {
            break;
        }
    }
    if (matrix_needed_ > 0) {
        throw std::runtime_error("CONTROL: incomplete MATRIX at EOF");
    }
    return results;
}

std::vector<ControlGoResult> ControlSession::run_script(
    const std::filesystem::path& script,
    const std::optional<std::filesystem::path>& image_override
) {
    std::ifstream in(script);
    if (!in) {
        throw std::runtime_error("CONTROL: cannot open " + script.string());
    }
    return run_stream(in, image_override);
}

}  // namespace mosflm::control
