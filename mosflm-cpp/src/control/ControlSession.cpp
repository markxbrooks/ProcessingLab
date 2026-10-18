#include "mosflm/control/ControlSession.hpp"

#include "mosflm/crystal/SpaceGroupInfo.hpp"
#include "mosflm/image/ImageIO.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace mosflm::control {
namespace {

AutoIndexRunner& autoindex_runner_slot() {
    static AutoIndexRunner runner;
    return runner;
}

std::string strip_comment(std::string_view line) {
    std::string out;
    out.reserve(line.size());
    for (char c : line) {
        if (c == '!' || c == '#') {
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
        "REFIX", "FINDSPOTS", "PICK", "XDISPLAY", "TCPIP", "SOCKET",
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
        for (std::size_t i = 1; i < args.size(); ++i) {
            const std::string u = upper(args[i]);
            if (u == "TO" && i + 1 < args.size()) {
                state_.process_last = require_int(args[i + 1], "PROCESS");
                ++i;
            } else if (u == "START" && i + 1 < args.size()) {
                state_.process_start_phi_deg =
                    require_double(args[i + 1], "PROCESS");
                ++i;
            } else if ((u == "ANGLE" || u == "OSCILLATION") &&
                       i + 1 < args.size()) {
                state_.process_angle_deg =
                    require_double(args[i + 1], "PROCESS");
                ++i;
            }
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
    if (key == "POSTREFINE" || key == "POSTREFL") {
        state_.enable_postrefinement = true;
        if (!args.empty() && (upper(args[0]) == "OFF" || args[0] == "0")) {
            state_.enable_postrefinement = false;
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
    if (key == "AUTOINDEX") {
        parse_autoindex_args(args);
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

std::filesystem::path ControlSession::resolve_image_path(
    const std::optional<std::filesystem::path>& image_override
) const {
    if (image_override && !image_override->empty()) {
        return *image_override;
    }
    if (!state_.image_path.empty()) {
        if (state_.image_path.is_absolute()) {
            return state_.image_path;
        }
        return state_.directory / state_.image_path;
    }
    if (!state_.image_template.empty()) {
        const std::string name =
            expand_template(state_.image_template, state_.process_first);
        return state_.directory / name;
    }
    throw std::runtime_error(
        "CONTROL: no IMAGE/TEMPLATE and no image override for PROCESS"
    );
}

diffraction::ProcessRequest ControlSession::make_process_request(
    const image::DiffractionImage& image
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
    // Match ProcessPipeline::known_indexed_insulin detector frame.
    det.yscale = state_.yscale;
    det.tilt_mrad = state_.tilt_mrad;
    det.twist_mrad = state_.twist_mrad;
    det.use_rotmat2 = state_.use_rotmat2;

    diffraction::ScanRange scan;
    scan.phi_start_deg = state_.process_start_phi_deg;
    scan.phi_end_deg =
        state_.process_start_phi_deg + state_.process_angle_deg;

    diffraction::PredictionSettings prediction;
    prediction.mosaic_spread_deg = state_.mosaic_deg;
    prediction.divh_deg = state_.has_divergence ? state_.divh_deg : 0.0;
    prediction.divv_deg = state_.has_divergence ? state_.divv_deg : 0.0;
    prediction.dmin_A = state_.dmin_A;
    prediction.dmax_A = state_.dmax_A;
    prediction.phi_step_deg = state_.phi_step_deg;
    prediction.apply_systematic_absences = true;
    // REEK-style gates tuned for insulin PROCESS (~1680 gen target).
    prediction.exclude_ice_rings = false;
    prediction.exclude_cusp = true;
    prediction.cusp_fraction = 0.035;
    prediction.max_phi_width_deg = 5.0;
    prediction.apply_reek_detector_margin = true;
    prediction.detector_margin_px = 24.0;
    prediction.widen_ewald_with_mosaic = true;
    prediction.mosaic_ewald_scale = 0.55;

    diffraction::IntegrationSettings integration;
    integration.box_half_px = 3;
    integration.bg_frame_px = 2;
    integration.apply_lp = true;
    integration.mosflm_detector_axes = true;
    integration.image_beam_x_mm =
        state_.has_beam ? state_.beam_x_mm : image.metadata.beam_x_mm;
    integration.image_beam_y_mm =
        state_.has_beam ? state_.beam_y_mm : image.metadata.beam_y_mm;

    diffraction::MeasurementBoxSettings measurement;
    measurement.separation_x_mm = state_.separation_x_mm;
    measurement.separation_y_mm = state_.separation_y_mm;
    measurement.bg_rim_px = 2;
    measurement.bg_reject_nsigma = 3.0;
    measurement.mark_overlaps = true;
    measurement.integrate_overlaps = true;

    diffraction::StandardProfileSettings profile;
    profile.box = measurement;
    profile.learn_min_isigi = 5.0;
    profile.min_learners = 8;
    profile.fallback_sigma_px = 1.5;
    profile.prefer_profile = true;
    profile.n_bins_x = state_.nmasks_x;
    profile.n_bins_y = state_.nmasks_y;

    refinement::PostRefinementSettings postref;

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
        state_.enable_postrefinement,
        postref,
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

    if (state_.autoindex_pending) {
        if (!autoindex_runner_slot()) {
            throw std::runtime_error(
                "CONTROL AUTOINDEX: DPS indexer not linked into this build"
            );
        }
        autoindex_runner_slot()(state_);
        state_.autoindex_pending = false;
        state_.done_run = true;
        out.ran_autoindex = true;
        out.message = "AUTOINDEX complete";
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
        out.message = "GO with no pending PROCESS/AUTOINDEX/STRATEGY";
        return out;
    }

    const auto img_path = resolve_image_path(image_override);
    out.image_used = img_path;
    // HDF5/NeXus masters: PROCESS first image is the 1-based frame index.
    const auto image =
        image::read_diffraction_image(img_path, state_.process_first);
    const auto req = make_process_request(image);
    out.crystal = req.crystal;
    out.process = diffraction::ProcessPipeline{}.run(req);

    // Multi-image NADD scaffolding: re-run adjacent frames when nadd > 1 and
    // TEMPLATE is set (joint postref still single-image refine for now).
    if (state_.nadd > 1 && !state_.image_template.empty() &&
        state_.enable_postrefinement) {
        // Observations from the primary process remain; flag in message.
        out.message = "PROCESS complete (NADD=" + std::to_string(state_.nadd) +
                      " reserved for joint refine)";
    } else {
        out.message = "PROCESS complete";
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
