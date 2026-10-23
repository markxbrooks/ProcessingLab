#include "mosflm/io/RunRecord.hpp"

#include "mosflm/crystal/SpaceGroupInfo.hpp"
#include "mosflm/diffraction/DataQualityStats.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <system_error>

namespace mosflm::io {
namespace {

/// RFC 3339 local time with offset, e.g. ``2026-10-23T09:36:59-07:00``.
std::string now_rfc3339() {
    const std::time_t t = std::time(nullptr);
    std::tm tm{};
    localtime_r(&t, &tm);
    char buf[40];
    std::strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%S%z", &tm);
    std::string s = buf;
    if (s.size() >= 5) {
        s.insert(s.size() - 2, ":");  // %z gives -0700
    }
    return s;
}

/// Wall-clock seconds since ``start`` (steady clock), to the millisecond.
double seconds_since(std::chrono::steady_clock::time_point start) {
    const auto d = std::chrono::steady_clock::now() - start;
    return std::round(std::chrono::duration<double>(d).count() * 1000.0) / 1000.0;
}

/// A stage object: every key present, optional ones null until set.
Json stage_json(const std::string& name) {
    return Json(Json::Object{{"name", name},           {"status", "pending"},
                             {"started_at", Json()},   {"finished_at", Json()},
                             {"duration_s", Json()},   {"message", Json()},
                             {"results", Json()},      {"error", Json()}});
}

Json space_group_json(int number) {
    const auto* sg = crystal::SpaceGroupInfo::find_by_number(number);
    return Json(Json::Object{{"number", number}, {"symbol", sg != nullptr ? Json(sg->hm) : Json()}});
}

std::string upper_first_token(const std::string& line) {
    std::string s = line.substr(0, line.find('!'));
    std::istringstream in(s);
    std::string tok;
    in >> tok;
    std::transform(tok.begin(), tok.end(), tok.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
    return tok;
}

Json cell_json(const UnitCell& c) {
    Json j = Json::object();
    j["a"] = c.a;
    j["b"] = c.b;
    j["c"] = c.c;
    j["alpha"] = c.alpha_deg;
    j["beta"] = c.beta_deg;
    j["gamma"] = c.gamma_deg;
    j["length_unit"] = "angstrom";
    j["angle_unit"] = "degree";
    return j;
}

/// A ``float`` setting as the decimal it was typed as (no 0.6999999881).
double f2d(float v) {
    return std::round(static_cast<double>(v) * 1e6) / 1e6;
}

Json shell_json(const diffraction::ResolutionShellStats& s) {
    // XDS's −99.9 "missing" markers become null.
    const auto opt = [](double v) { return v <= -99.0 ? Json() : Json(v); };
    Json j = Json::object();
    j["d_min_angstrom"] = s.d_limit_A;
    j["observed"] = s.observed;
    j["unique"] = s.unique;
    j["possible"] = s.possible;
    j["completeness_pct"] = s.completeness_pct;
    j["r_obs_pct"] = opt(s.r_obs_pct);
    j["r_meas_pct"] = opt(s.rmeas_pct);
    j["mean_i_over_sigma"] = opt(s.mean_isigma);
    j["cc_half_pct"] = s.compared > 0 ? Json(s.cc_half_pct) : Json();
    j["compared"] = s.compared;
    return j;
}

/// A resolved ``input.json`` experiment value, or nullptr.
const Json* input_value(const Json& input, const char* key) {
    const Json* e = input.find("experiment");
    const Json* f = e != nullptr ? e->find(key) : nullptr;
    const Json* v = f != nullptr ? f->find("value") : nullptr;
    return v != nullptr && !v->is_null() ? v : nullptr;
}

/// The run's final experiment model: refined values from ``st`` where
/// CONTROL holds them, the resolved inputs otherwise (CONTROL defaults are
/// not authoritative). The orientation carries the cell its matrix implies.
Json final_experiment(const Json& input, const control::ControlState& st, int first, int last) {
    const auto pick = [&](bool held, double state_value, const char* key, std::size_t k) {
        if (held) {
            return state_value;
        }
        const Json* v = input_value(input, key);
        if (v == nullptr) {
            return state_value;
        }
        return v->is_array() ? v->items()[k].number() : v->number();
    };
    const double lam = pick(st.has_wavelength, st.wavelength_A, "wavelength", 0);
    Json x = Json::object();
    x["crystal"] = Json(Json::Object{{"unit_cell", cell_json(st.cell)},
                                     {"space_group", st.has_symmetry ? space_group_json(st.symmetry_number) : Json()}});
    x["wavelength_angstrom"] = lam;
    Json d = Json::object();
    d["distance_mm"] = pick(st.has_distance, st.distance_mm, "detector_distance", 0);
    d["beam_centre_mm"] = Json(Json::Object{{"x", pick(st.has_beam, st.beam_x_mm, "beam_centre", 0)},
                                            {"y", pick(st.has_beam, st.beam_y_mm, "beam_centre", 1)}});
    d["pixel_size_mm"] = Json(Json::Object{{"x", pick(st.has_pixel, st.pixel_x_mm, "pixel_size", 0)},
                                           {"y", pick(st.has_pixel, st.pixel_y_mm, "pixel_size", 1)}});
    d["tilt_twist_mrad"] = Json(Json::Object{{"tilt", st.tilt_mrad}, {"twist", st.twist_mrad}});
    if (const Json* e = input.find("experiment"); e != nullptr && e->find("detector_size_px") != nullptr) {
        d["size_px"] = *e->find("detector_size_px");
    }
    x["detector"] = std::move(d);
    Json scan = Json::object();
    scan["frames"] = Json(Json::Object{{"first", first}, {"last", last}});
    const Json* osc = input_value(input, "oscillation_width");
    scan["oscillation_deg"] = osc != nullptr ? *osc : Json();
    x["scan"] = std::move(scan);
    x["mosaicity_deg"] = st.mosaic_deg;
    if (st.has_matrix) {
        Json rows = Json::array();
        for (int i = 0; i < 3; ++i) {
            rows.push_back(Json(Json::Array{Json(st.amat(i, 0)), Json(st.amat(i, 1)), Json(st.amat(i, 2))}));
        }
        const UnitCell implied = OrientationMatrix::from_mosflm_orientation_matrix(st.amat, lam).cell();
        Json o = Json::object();
        o["convention"] = "mosflm_amat";
        o["matrix"] = std::move(rows);
        o["wavelength_angstrom"] = lam;
        o["invertx"] = st.spt_invertx;
        o["unit_cell"] = cell_json(implied);
        x["orientation"] = std::move(o);
    } else {
        x["orientation"] = Json();
    }
    x["conventions"] = Json(Json::Object{
        {"orientation", "mosflm_amat: MOSFLM [A] = λ·U·B; rows are the MOSFLM lab axes X (along the beam), Y, "
                        "Z (spindle); columns are a*, b*, c*; scaled by orientation.wavelength_angstrom. "
                        "invertx is the detector handedness it was derived in (CONTROL INVERTX). "
                        "orientation.unit_cell is the cell the matrix implies."},
        {"beam_centre", "x, y in mm along ProcessingLab's detector axes: header pixel coordinates × pixel size "
                        "(CBF/HDF5: image fast, slow; ADSC SMV: stored transposed, as MOSFLM reads it). "
                        "XDS.INP ORGX×QX, ORGY×QY and CONTROL BEAM use the same frame."}});
    return x;
}

std::string fmt(double v, int prec) {
    std::ostringstream o;
    o << std::fixed << std::setprecision(prec) << v;
    return o.str();
}

/// Markdown value of an experiment field's ``value``.
std::string md_value(const Json& v) {
    if (v.is_null()) {
        return "—";
    }
    if (v.is_number()) {
        char buf[32];
        std::snprintf(buf, sizeof buf, "%.6g", v.number());
        return buf;
    }
    if (v.is_array()) {
        std::string s;
        for (const auto& e : v.items()) {
            s += (s.empty() ? "" : " ") + md_value(e);
        }
        return s;
    }
    if (v.is_object()) {
        if (const auto* f = v.find("first")) {
            return md_value(*f) + "–" + md_value(*v.find("last"));
        }
    }
    return v.string();
}

std::string md_num(const Json* v, int prec, const std::string& suffix = "") {
    return v == nullptr || v->is_null() ? "—" : fmt(v->number(), prec) + suffix;
}

}  // namespace

std::string pending_stage(const control::ControlState& st) {
    if (st.findspots_pending) {
        return "findspots";
    }
    if (st.autoindex_pending) {
        return "autoindex";
    }
    if (st.process_pending) {
        return "process";
    }
    if (st.strategy_pending) {
        return "strategy";
    }
    return "go";
}

std::vector<std::string> planned_stages(const std::string& control_text) {
    // GO runs the first pending step in this order (ControlSession::go).
    static const std::vector<std::pair<std::string, std::string>> order{
        {"FINDSPOTS", "findspots"}, {"AUTOINDEX", "autoindex"},
        {"PROCESS", "process"}, {"STRATEGY", "strategy"}};
    std::vector<std::string> pending;
    std::vector<std::string> stages;
    std::istringstream in(control_text);
    std::string line;
    while (std::getline(in, line)) {
        const std::string key = upper_first_token(line);
        if (key == "EXIT" || key == "QUIT" || key == "END") {
            break;
        }
        if (key == "GO" || key == "RUN") {
            for (const auto& [kw, name] : order) {
                const auto it = std::find(pending.begin(), pending.end(), name);
                if (it != pending.end()) {
                    stages.push_back(name);
                    pending.erase(it);
                    break;
                }
            }
            continue;
        }
        for (const auto& [kw, name] : order) {
            if (key == kw && std::find(pending.begin(), pending.end(), name) == pending.end()) {
                pending.push_back(name);
            }
        }
    }
    return stages;
}

Json settings_json(const control::ControlState& st) {
    Json s = Json::object();
    const auto& fs = st.findspots;
    Json f = Json::object();
    f["threshold_sigma"] = f2d(fs.threshold_sigma);
    f["max_spots"] = fs.max_spots;
    f["frames"] = fs.multiframe
                      ? Json(Json::Array{Json(fs.range_first), Json(fs.range_last)})
                      : Json(Json::Array{Json(fs.image_number), Json(fs.image_number)});
    f["merge_xy_mm"] = fs.merge_xy_mm;
    f["spot_file"] = fs.spot_file.string();
    s["findspots"] = std::move(f);
    const auto& ai = st.autoindex;
    Json a = Json::object();
    a["threshold_i_over_sigma"] = f2d(ai.min_isigi);
    a["dmax_angstrom"] = st.dmax_A;
    a["methods"] = Json(Json::Array{});
    if (ai.use_dps) a["methods"].push_back("dps");
    if (ai.use_fft3d) a["methods"].push_back("fft3d");
    if (ai.use_diff_vector) a["methods"].push_back("difference_vectors");
    a["sepmin"] = f2d(ai.sepmin);
    a["cluster_radius"] = f2d(ai.cluster_radius);
    a["geometry_refinement"] = ai.use_idxref_refine;
    a["min_indexed_fraction"] = f2d(ai.min_indexed_frac);
    s["autoindex"] = std::move(a);
    Json i = Json::object();
    i["frames"] = Json(Json::Array{Json(st.process_first), Json(st.process_last)});
    i["dmin_angstrom"] = st.dmin_A > 0.0 ? Json(st.dmin_A) : Json("detector edge");
    i["peak_box"] = st.peak_half_px >= 1 ? Json(2 * st.peak_half_px + 1) : Json("from spot size");
    i["profile"] = st.profile_3d ? "3d" : "2d";
    i["profile_regions"] = Json(Json::Array{Json(st.nmasks_x), Json(st.nmasks_y)});
    i["rocking_model"] = st.rocking_model;
    i["mosaicity_deg"] = st.mosaic_from_user ? Json(st.mosaic_deg) : Json("estimated");
    i["partiality_scaling"] = st.scale_by_partiality;
    s["integration"] = std::move(i);
    Json r = Json::object();
    r["period_frames"] = st.postrefine_period;
    Json what = Json::array();
    if (st.refine_orientation) what.push_back("orientation");
    if (st.refine_cell) what.push_back("cell");
    if (st.refine_distance) what.push_back("distance");
    if (st.refine_beam) what.push_back("beam");
    if (st.refine_tilt_twist) what.push_back("tilt_twist");
    r["parameters"] = std::move(what);
    s["refinement"] = std::move(r);
    s["threads"] = st.nproc;
    return s;
}

Json input_record(
    const ExperimentModel& m, const control::ControlState& st, const RunContext& ctx
) {
    Json j = Json::object();
    j["schema_version"] = kRunRecordSchemaVersion;
    j["program"] = Json(Json::Object{{"name", "ProcessingLab"}, {"version", "0.1.0"}});
    j["mode"] = ctx.mode;
    j["run_mode"] = run_mode_name(ctx.run_mode);
    j["parameter_source"] = run_mode_parameter_source(ctx.run_mode);
    j["created_at"] = now_rfc3339();
    Json cl = Json::array();
    for (const auto& a : ctx.command_line) {
        cl.push_back(a);
    }
    j["command_line"] = std::move(cl);
    Json inputs = Json::object();
    Json img = Json::object();
    img["path"] = m.image.value.string();
    img["source"] = source_id(m.image.source);
    img["format"] = m.hdf5 ? "hdf5" : "frame_files";
    if (!m.hdf5) {
        img["template"] = (m.directory / m.image_template).string();
    }
    img["frames_available"] = m.frames_available;
    inputs["images"] = Json(Json::Array{img});
    Json imp = Json::object();
    imp["xds_inp"] = ctx.xds_inp ? Json(ctx.xds_inp->string()) : Json();
    imp["xinfo"] = ctx.xinfo ? Json(ctx.xinfo->string()) : Json();
    imp["xinfo_sweep"] = ctx.sweep ? Json(*ctx.sweep) : Json();
    imp["xparm"] = ctx.xparm ? Json(ctx.xparm->string()) : Json();
    imp["previous_run"] = ctx.previous_run ? Json(ctx.previous_run->string()) : Json();
    inputs["imports"] = std::move(imp);
    j["inputs"] = std::move(inputs);
    j["experiment"] = m.to_json();
    j["settings"] = settings_json(st);
    Json w = Json::array();
    for (const auto& s : m.warnings) {
        w.push_back(s);
    }
    j["warnings"] = std::move(w);
    j["control"] = ctx.control_text;
    return j;
}

ProcessRecord::ProcessRecord(
    std::filesystem::path dir, Json input, const std::vector<std::string>& stages
)
    : dir_(std::move(dir)), input_(std::move(input)) {
    doc_["schema_version"] = kRunRecordSchemaVersion;
    doc_["status"] = "running";
    doc_["started_at"] = now_rfc3339();
    doc_["finished_at"] = Json();
    doc_["duration_s"] = Json();
    doc_["run_directory"] = std::filesystem::absolute(dir_).lexically_normal().string();
    // Operational completion is not scientific validation; no checks run yet.
    doc_["validation"] = Json(Json::Object{{"status", "not_assessed"}, {"checks", Json::array()}});
    Json st = Json::array();
    for (const auto& name : stages) {
        st.push_back(stage_json(name));
    }
    doc_["stages"] = std::move(st);
    doc_["experiment"] = Json();
    doc_["indexing"] = Json();
    doc_["refinement"] = Json();
    doc_["integration"] = Json();
    doc_["quality"] = Json();
    doc_["warnings"] = Json::array();
    for (const auto& w : input_["warnings"].items()) {
        doc_["warnings"].push_back(w);
    }
    doc_["error"] = Json();
    doc_["outputs"] = Json::array();
    write_file_atomic(dir_ / "input.json", input_.dump());
    write();
}

Json& ProcessRecord::stage(const std::string& name) {
    // The first stage of that name not yet finished; unplanned ones are added.
    auto& arr = doc_["stages"].elements();
    for (auto& s : arr) {
        const std::string status = s.find("status") ? s.find("status")->string() : "";
        if (s.find("name") && s.find("name")->string() == name &&
            (status == "pending" || status == "running")) {
            return s;
        }
    }
    arr.push_back(stage_json(name));
    return arr.back();
}

void ProcessRecord::begin_stage(const control::ControlState& state) {
    running_ = pending_stage(state);
    Json& s = stage(running_);
    s["status"] = "running";
    s["started_at"] = now_rfc3339();
    stage_start_ = std::chrono::steady_clock::now();
    write();
}

void ProcessRecord::end_stage(const control::ControlGoResult& r, const control::ControlState& st) {
    Json& s = stage(running_);
    s["status"] = "completed";
    s["finished_at"] = now_rfc3339();
    s["duration_s"] = seconds_since(stage_start_);
    if (!r.message.empty()) {
        s["message"] = r.message;
    }
    if (r.ran_findspots) {
        Json res = Json::object();
        res["spots"] = st.findspots_n_spots;
        s["results"] = std::move(res);
        add_output(control::output_path(st, st.findspots.spot_file), "spot_list");
    }
    if (r.ran_autoindex) {
        // Stage summary; the full solution is in ``indexing``.
        s["results"] = Json(Json::Object{{"lattice", st.autoindex_system},
                                         {"indexed_fraction", st.autoindex_indexed_frac},
                                         {"rmsd_hkl", st.autoindex_rmsd_hkl}});
        Json ix = Json::object();
        ix["lattice"] = st.autoindex_system;
        ix["lattice_character"] = st.autoindex_character;
        ix["space_group"] = st.has_symmetry ? space_group_json(st.symmetry_number) : Json();
        ix["unit_cell"] = st.has_cell ? cell_json(st.cell) : Json();
        Json mat = Json::array();
        for (int i = 0; i < 3; ++i) {
            mat.push_back(Json(Json::Array{Json(st.amat(i, 0)), Json(st.amat(i, 1)), Json(st.amat(i, 2))}));
        }
        ix["orientation_matrix_mosflm"] = st.has_matrix ? std::move(mat) : Json();
        ix["spots_used"] = st.autoindex_n_spots;
        ix["indexed_fraction"] = st.autoindex_indexed_frac;
        ix["rmsd_hkl"] = st.autoindex_rmsd_hkl;
        doc_["indexing"] = Json(Json::Object{{"selected_solution", std::move(ix)}});
        add_output(control::output_path(st, "XPARM.XDS"), "xparm");
        add_output(control::output_path(st, "index.log"), "index_summary");
    }
    if (r.ran_process) {
        // The final model, after the last periodic refinement of PROCESS.
        doc_["experiment"] = final_experiment(input_, st, r.process_first, r.process_last);
        // What PROCESS did to each parameter (the policy it ran with).
        Json ref = Json::object();
        ref["stage"] = "process";
        Json params = Json::object();
        params["orientation"] = st.refine_orientation ? "refined" : "fixed";
        params["unit_cell"] = st.refine_cell ? "refined" : "fixed";
        params["detector_distance"] = st.refine_distance ? "refined" : "fixed";
        params["beam_centre"] = st.refine_beam ? "refined" : "fixed";
        params["detector_tilt_twist"] = st.refine_tilt_twist ? "refined" : "fixed";
        params["mosaicity"] = st.mosaic_from_user ? "fixed" : "estimated";
        ref["parameters"] = std::move(params);
        ref["unit_cell_constraint"] = st.has_symmetry ? "space_group_metric" : "none";
        doc_["refinement"] = std::move(ref);
        const auto& p = r.process;
        Json in = Json::object();
        in["frames"] = Json(Json::Object{{"first", r.process_first}, {"last", r.process_last}});
        const int half = st.peak_half_px >= 1 ? st.peak_half_px : st.spot_peak_half_px;
        in["peak_box_px"] = half >= 1 ? Json(2 * half + 1) : Json();
        // Populations: see docs/RUN_RECORDS.md. The first four count
        // predicted observations (per frame, before partials are summed);
        // integrated_reflections counts the rows written to the MTZ.
        Json counts = Json::object();
        counts["predicted_observations"] = p.predictions.size();
        counts["profile_fitted_observations"] = p.n_profile_ok;
        counts["overlapped_observations"] = p.n_overlapped;
        counts["omitted"] = p.n_omitted;
        counts["integrated_reflections"] = p.mosflm_reflections.size();
        in["counts"] = std::move(counts);
        doc_["integration"] = std::move(in);
        s["results"] = Json(Json::Object{{"integrated_reflections", p.mosflm_reflections.size()}});
        if (!p.mosflm_reflections.empty() && r.crystal.spacegroup != nullptr) {
            const auto q = diffraction::compute_data_quality(
                r.crystal, p.mosflm_reflections, st.dmax_A, st.dmin_A
            );
            Json qual = Json::object();
            qual["scaling_status"] = "unscaled";
            qual["outlier_rejection_status"] = "not_performed";
            qual["space_group"] = space_group_json(r.crystal.spacegroup->number);
            qual["resolution_range_angstrom"] = Json(Json::Object{{"low", q.dmax_A}, {"high", q.dmin_A}});
            qual["min_i_over_sigma"] = q.min_isigma;
            qual["shell_definition"] = Json(Json::Object{
                {"method", "equal_width_in_inverse_d"},
                {"count", q.shells.size()},
                {"d_min_angstrom", "high-resolution edge of each shell"}});
            Json shells = Json::array();
            for (const auto& sh : q.shells) {
                shells.push_back(shell_json(sh));
            }
            qual["shells"] = std::move(shells);
            qual["total"] = shell_json(q.total);
            // One population, the integrated reflections: these sum to it.
            Json qc = Json::object();
            qc["integrated_reflections"] = p.mosflm_reflections.size();
            qc["rejected_misfit_observations"] = q.n_rejected_misfits;
            qc["invalid_sigma_observations"] = q.n_invalid_sigma;
            qc["systematic_absence_observations"] = q.n_systematic_absent;
            qc["below_min_i_over_sigma_observations"] = q.n_below_min_isigma;
            qc["outside_resolution_observations"] = q.n_outside_resolution;
            qc["accepted_observations"] = q.n_accepted_obs;
            qc["unique_reflections"] = q.n_unique_accepted;
            qual["counts"] = std::move(qc);
            doc_["quality"] = std::move(qual);
        }
        if (st.has_xparm_out) {
            add_output(st.xparm_out, "gxparm");
        }
    }
    running_.clear();
    write();
}

void ProcessRecord::add_output(const std::filesystem::path& p, const std::string& kind) {
    std::error_code ec;
    if (!std::filesystem::exists(p, ec)) {
        return;
    }
    // Relative to the run directory when inside it, so a run can be moved.
    std::filesystem::path rel = p;
    const auto r = std::filesystem::absolute(p).lexically_normal().lexically_relative(
        std::filesystem::absolute(dir_).lexically_normal());
    if (!r.empty() && *r.begin() != "..") {
        rel = r;
    }
    for (const auto& o : doc_["outputs"].items()) {
        if (o.find("path") && o.find("path")->string() == rel.string()) {
            return;
        }
    }
    doc_["outputs"].push_back(Json(Json::Object{{"kind", kind}, {"path", rel.string()}}));
    write();
}

void ProcessRecord::add_warning(const std::string& w) {
    doc_["warnings"].push_back(w);
    write();
}

void ProcessRecord::fail(const std::string& error) {
    if (!running_.empty()) {
        Json& s = stage(running_);
        s["status"] = "failed";
        s["finished_at"] = now_rfc3339();
        s["duration_s"] = seconds_since(stage_start_);
        s["error"] = error;
    }
    doc_["status"] = "failed";
    doc_["error"] = error;
    doc_["finished_at"] = now_rfc3339();
    doc_["duration_s"] = seconds_since(run_start_);
    write();
}

void ProcessRecord::complete() {
    // Planned stages the CONTROL never reached stay visible as skipped.
    auto& arr = doc_["stages"].elements();
    for (auto& s : arr) {
        if (s.find("status") && s.find("status")->string() == "pending") {
            s["status"] = "skipped";
        }
    }
    doc_["status"] = "completed";
    doc_["finished_at"] = now_rfc3339();
    doc_["duration_s"] = seconds_since(run_start_);
    write();
}

void ProcessRecord::write() {
    write_file_atomic(dir_ / "process.json", doc_.dump());
    write_file_atomic(dir_ / "process.md", render_markdown(input_, doc_, dir_));
}

std::string render_markdown(const Json& in, const Json& pr, const std::filesystem::path& dir) {
    std::ostringstream md;
    const auto str = [](const Json* j) { return j == nullptr ? std::string() : j->string(); };
    md << "# ProcessingLab report\n\n";
    md << "**Status:** " << str(pr.find("status")) << "  \n";
    md << "**Started:** " << str(pr.find("started_at"));
    if (const auto* f = pr.find("finished_at"); f != nullptr && !f->is_null()) {
        md << " · **Finished:** " << f->string();
    }
    if (const auto* d = pr.find("duration_s"); d != nullptr && !d->is_null()) {
        md << " (" << fmt(d->number(), 1) << " s)";
    }
    md << "  \n**Mode:** `ProcessingLab " << (str(in.find("mode")) == "run" ? "run" : "<script>") << "`";
    if (const auto* rm = in.find("run_mode")) {
        md << " · **Run mode:** `" << rm->string() << "` — parameters from "
           << str(in.find("parameter_source"));
    }
    md << "\n\n";
    if (const auto* e = pr.find("error"); e != nullptr && !e->is_null()) {
        md << "> **Error:** " << e->string() << "\n\n";
    }

    md << "## Input\n\n";
    if (const auto* inputs = in.find("inputs")) {
        for (const auto& img : (*inputs).find("images")->items()) {
            md << "- Images: `" << str(img.find("path")) << "`";
            if (const auto* t = img.find("template")) {
                md << " (template `" << t->string() << "`)";
            }
            const auto n = img.find("frames_available")->number();
            md << ", " << md_value(*img.find("frames_available")) << (n == 1 ? " frame" : " frames")
               << " on disk\n";
        }
        for (const auto& [k, v] : inputs->find("imports")->members()) {
            if (!v.is_null()) {
                md << "- " << k << ": `" << v.string() << "`\n";
            }
        }
    }
    md << "\n| Parameter | Value | Source | Overrides |\n| --- | --- | --- | --- |\n";
    static const std::vector<std::pair<std::string, std::string>> labels{
        {"wavelength", "Wavelength (Å)"}, {"detector_distance", "Distance (mm)"},
        {"beam_centre", "Beam centre (mm)"}, {"pixel_size", "Pixel (mm)"},
        {"oscillation_width", "Oscillation (°)"}, {"unit_cell", "Unit cell"},
        {"space_group_number", "Space group"}, {"frame_range", "Frames"}};
    if (const auto* exp = in.find("experiment")) {
        for (const auto& [key, label] : labels) {
            const auto* f = exp->find(key);
            if (f == nullptr || f->find("value")->is_null()) {
                continue;
            }
            std::string ov;
            const std::string now = md_value(*f->find("value"));
            for (const auto& o : f->find("overrides")->items()) {
                const std::string was = md_value(*o.find("value"));
                ov += (ov.empty() ? "" : "; ") +
                      (was == now ? "same in " + str(o.find("source")) : was + " (" + str(o.find("source")) + ")");
            }
            md << "| " << label << " | " << md_value(*f->find("value")) << " | "
               << str(f->find("source")) << " | " << (ov.empty() ? "—" : ov) << " |\n";
        }
    }

    md << "\n## Stages\n\n| Stage | Status | Detail |\n| --- | --- | --- |\n";
    for (const auto& s : pr.find("stages")->items()) {
        std::string detail = str(s.find("message"));
        if (const auto* e = s.find("error"); e != nullptr && !e->is_null()) {
            detail = e->string();
        }
        if (const auto* r = s.find("results"); r != nullptr && r->find("spots")) {
            detail = md_value(*r->find("spots")) + " spots";
        }
        md << "| " << str(s.find("name")) << " | " << str(s.find("status")) << " | "
           << (detail.empty() ? "—" : detail) << " |\n";
    }

    if (const auto* ix = pr.find("indexing"); ix != nullptr && !ix->is_null()) {
        const Json& sol = *ix->find("selected_solution");
        md << "\n## Indexing\n\n";
        md << "- Lattice: " << str(sol.find("lattice"));
        if (const auto* sg = sol.find("space_group"); sg && !sg->is_null()) {
            md << ", space group " << md_value(*sg->find("number"));
        }
        md << "\n";
        if (const auto* c = sol.find("unit_cell"); c && !c->is_null()) {
            md << "- Cell: " << md_num(c->find("a"), 2) << " " << md_num(c->find("b"), 2) << " "
               << md_num(c->find("c"), 2) << " Å, " << md_num(c->find("alpha"), 2) << " "
               << md_num(c->find("beta"), 2) << " " << md_num(c->find("gamma"), 2) << "°\n";
        }
        md << "- Indexed " << fmt(100.0 * sol.find("indexed_fraction")->number(), 0) << "% of "
           << md_value(*sol.find("spots_used")) << " spots, RMSD " << md_num(sol.find("rmsd_hkl"), 3)
           << " (hkl)\n";
    }

    if (const auto* x = pr.find("experiment"); x != nullptr && !x->is_null()) {
        const Json& c = *x->find("crystal")->find("unit_cell");
        const Json& det = *x->find("detector");
        const Json& beam = *det.find("beam_centre_mm");
        md << "\n## Refined geometry\n\n";
        md << "- Cell: " << md_num(c.find("a"), 3) << " " << md_num(c.find("b"), 3) << " "
           << md_num(c.find("c"), 3) << " Å, " << md_num(c.find("alpha"), 2) << " "
           << md_num(c.find("beta"), 2) << " " << md_num(c.find("gamma"), 2) << "°\n";
        md << "- Distance: " << md_num(det.find("distance_mm"), 3) << " mm; beam "
           << md_num(beam.find("x"), 3) << " " << md_num(beam.find("y"), 3) << " mm\n";
        md << "- Mosaicity: " << md_num(x->find("mosaicity_deg"), 3) << "°\n";
        if (const auto* rf = pr.find("refinement"); rf != nullptr && !rf->is_null()) {
            std::string refined;
            std::string fixed;
            for (const auto& [k, v] : rf->find("parameters")->members()) {
                auto& list = v.string() == "fixed" ? fixed : refined;
                list += (list.empty() ? "" : ", ") + k;
            }
            md << "- Refined: " << (refined.empty() ? "nothing" : refined) << "; fixed: "
               << (fixed.empty() ? "nothing" : fixed) << "\n";
        }
    }

    if (const auto* ig = pr.find("integration"); ig != nullptr && !ig->is_null()) {
        md << "\n## Integration\n\n";
        md << "- Frames " << md_value(*ig->find("frames")) << ", peak box "
           << md_value(*ig->find("peak_box_px")) << " px\n";
        const Json& c = *ig->find("counts");
        md << "- " << md_value(*c.find("integrated_reflections")) << " reflections integrated from "
           << md_value(*c.find("predicted_observations")) << " predicted observations ("
           << md_value(*c.find("profile_fitted_observations")) << " profile-fitted, "
           << md_value(*c.find("overlapped_observations")) << " overlapped)\n";
    }

    if (const auto* q = pr.find("quality"); q != nullptr && !q->is_null()) {
        md << "\n## Data quality\n\n" << str(q->find("space_group")->find("symbol"))
           << "; unscaled: no batch scaling or outlier rejection.\n\n"
           << "| d_min (Å) | Observed | Unique | Completeness | R_meas | I/σ | CC½ |\n"
           << "| ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n";
        const auto row = [&](const Json& sh, const std::string& label) {
            md << "| " << label << " | " << md_value(*sh.find("observed")) << " | "
               << md_value(*sh.find("unique")) << " | " << md_num(sh.find("completeness_pct"), 1, "%")
               << " | " << md_num(sh.find("r_meas_pct"), 1, "%") << " | "
               << md_num(sh.find("mean_i_over_sigma"), 2) << " | " << md_num(sh.find("cc_half_pct"), 1)
               << " |\n";
        };
        for (const auto& sh : q->find("shells")->items()) {
            row(sh, fmt(sh.find("d_min_angstrom")->number(), 2));
        }
        row(*q->find("total"), "**total**");
        const Json& c = *q->find("counts");
        md << "\n" << md_value(*c.find("accepted_observations")) << " of "
           << md_value(*c.find("integrated_reflections")) << " integrated reflections accepted ("
           << md_value(*c.find("rejected_misfit_observations")) << " misfits, "
           << md_value(*c.find("systematic_absence_observations")) << " systematic absences, "
           << md_value(*c.find("below_min_i_over_sigma_observations")) << " below I/σ "
           << md_num(q->find("min_i_over_sigma"), 1) << ", "
           << md_value(*c.find("outside_resolution_observations")) << " outside the resolution range, "
           << md_value(*c.find("invalid_sigma_observations")) << " with invalid σ).\n";
    }

    if (const auto* w = pr.find("warnings"); w != nullptr && !w->items().empty()) {
        md << "\n## Warnings\n\n";
        for (const auto& s : w->items()) {
            md << "- " << s.string() << "\n";
        }
    }
    md << "\n## Outputs\n\n";
    for (const auto& o : pr.find("outputs")->items()) {
        std::filesystem::path path = str(o.find("path"));
        // Inside the run directory (xds/XDS.INP too): relative to it.
        if (!dir.empty()) {
            const auto rel = path.lexically_relative(dir);
            if (!rel.empty() && *rel.begin() != "..") {
                path = rel;
            }
        }
        md << "- " << str(o.find("kind")) << ": `" << path.string() << "`\n";
    }
    md << "\nStructured records: `input.json`, `process.json`. Full diagnostics: `process.log`.\n";
    return md.str();
}

namespace {

/// ``process.json`` ``experiment`` → the model ``--from`` carries forward,
/// refusing missing geometry, a cell the space group forbids, a matrix that
/// implies another cell, and a wavelength or pixel size that disagrees with
/// the run's ``input.json``.
PreviousModel read_final_model(const Json& x, const Json& input, const std::filesystem::path& pj) {
    const auto need = [&](const Json* parent, const char* key) -> const Json& {
        const Json* v = parent != nullptr ? parent->find(key) : nullptr;
        if (v == nullptr || v->is_null()) {
            throw std::runtime_error(pj.string() + ": experiment." + key + " is missing");
        }
        return *v;
    };
    const auto cell_of = [&](const Json& c) {
        return std::array<double, 6>{need(&c, "a").number(),     need(&c, "b").number(),
                                     need(&c, "c").number(),     need(&c, "alpha").number(),
                                     need(&c, "beta").number(),  need(&c, "gamma").number()};
    };
    PreviousModel m;
    const Json& crystal = need(&x, "crystal");
    m.cell = cell_of(need(&crystal, "unit_cell"));
    m.space_group = static_cast<int>(need(&need(&crystal, "space_group"), "number").number());
    check_cell_fits_spacegroup(m.cell, m.space_group);
    m.wavelength_A = need(&x, "wavelength_angstrom").number();
    const Json& det = need(&x, "detector");
    m.distance_mm = need(&det, "distance_mm").number();
    const Json& beam = need(&det, "beam_centre_mm");
    m.beam_mm = {need(&beam, "x").number(), need(&beam, "y").number()};
    const Json& pixel = need(&det, "pixel_size_mm");
    m.pixel_mm = {need(&pixel, "x").number(), need(&pixel, "y").number()};
    const Json& tt = need(&det, "tilt_twist_mrad");
    m.tilt_twist_mrad = {need(&tt, "tilt").number(), need(&tt, "twist").number()};
    m.mosaicity_deg = need(&x, "mosaicity_deg").number();

    const Json& o = need(&x, "orientation");
    if (need(&o, "convention").string() != "mosflm_amat") {
        throw std::runtime_error(pj.string() + ": orientation convention '" + o.find("convention")->string() +
                                 "' is not mosflm_amat");
    }
    const auto& rows = need(&o, "matrix").items();
    if (rows.size() != 3) {
        throw std::runtime_error(pj.string() + ": orientation.matrix is not 3×3");
    }
    Eigen::Matrix3d a;
    for (int r = 0; r < 3; ++r) {
        if (rows[static_cast<std::size_t>(r)].items().size() != 3) {
            throw std::runtime_error(pj.string() + ": orientation.matrix is not 3×3");
        }
        for (int k = 0; k < 3; ++k) {
            a(r, k) = rows[static_cast<std::size_t>(r)].items()[static_cast<std::size_t>(k)].number();
            m.matrix[static_cast<std::size_t>(3 * r + k)] = a(r, k);
        }
    }
    m.matrix_wavelength_A = need(&o, "wavelength_angstrom").number();
    m.invertx = need(&o, "invertx").boolean();
    // The matrix must belong to the cell it is carried with: recompute the
    // cell it implies rather than trusting the stored copy.
    const UnitCell implied = OrientationMatrix::from_mosflm_orientation_matrix(a, m.matrix_wavelength_A).cell();
    const std::array<double, 6> ic{implied.a, implied.b, implied.c,
                                   implied.alpha_deg, implied.beta_deg, implied.gamma_deg};
    for (std::size_t k = 0; k < 6; ++k) {
        const bool bad = k < 3 ? std::fabs(ic[k] - m.cell[k]) > 2e-3 * m.cell[k] : std::fabs(ic[k] - m.cell[k]) > 0.2;
        if (bad) {
            std::ostringstream s;
            s << std::fixed << std::setprecision(3) << pj.string() << ": the orientation matrix implies cell "
              << ic[0] << ' ' << ic[1] << ' ' << ic[2] << ' ' << ic[3] << ' ' << ic[4] << ' ' << ic[5]
              << ", not the run's " << m.cell[0] << ' ' << m.cell[1] << ' ' << m.cell[2] << ' ' << m.cell[3] << ' '
              << m.cell[4] << ' ' << m.cell[5];
            throw std::runtime_error(s.str());
        }
    }
    // The beam and detector are the run's own: they must match its inputs.
    const auto agrees = [&](const char* key, std::size_t k, double v, const char* what) {
        const Json* in = input_value(input, key);
        if (in == nullptr) {
            return;
        }
        const double want = in->is_array() ? in->items()[k].number() : in->number();
        if (std::fabs(v - want) > 1e-4 * std::fabs(want)) {
            throw std::runtime_error(pj.string() + ": final " + what + " " + fmt(v, 6) + " differs from " + fmt(want, 6) +
                                     " in input.json");
        }
    };
    agrees("wavelength", 0, m.wavelength_A, "wavelength");
    agrees("pixel_size", 0, m.pixel_mm.first, "pixel size x");
    agrees("pixel_size", 1, m.pixel_mm.second, "pixel size y");
    return m;
}

}  // namespace

PreviousRun read_previous_run(const std::filesystem::path& path) {
    namespace fs = std::filesystem;
    PreviousRun prev;
    prev.directory = fs::is_directory(path) ? path : path.parent_path();
    prev.directory = fs::absolute(prev.directory).lexically_normal();
    const fs::path pj = prev.directory / "process.json";
    if (!fs::exists(pj)) {
        std::string sweeps;
        std::error_code ec;
        for (const auto& e : fs::directory_iterator(prev.directory, ec)) {
            const std::string n = e.path().filename().string();
            if (e.is_directory() && n.rfind("sweep_", 0) == 0 && fs::exists(e.path() / "process.json")) {
                sweeps += (sweeps.empty() ? "" : ", ") + n;
            }
        }
        if (!sweeps.empty()) {
            throw std::runtime_error(prev.directory.string() + " holds several sweeps (" + sweeps +
                                     "); give one, e.g. --from " + (prev.directory / "sweep_1").string());
        }
        throw std::runtime_error("no process.json in " + prev.directory.string());
    }
    const Json process = read_json_file(pj);
    // Version 1 (before RFC 3339 / unit fields) and 2 are both readable.
    if (const auto* v = process.find("schema_version");
        v != nullptr && v->number() > kRunRecordSchemaVersion) {
        throw std::runtime_error(pj.string() + ": schema_version " + md_value(*v) +
                                 " is newer than this ProcessingLab reads (" +
                                 std::to_string(kRunRecordSchemaVersion) + ")");
    }
    const auto* status = process.find("status");
    if (status == nullptr || status->string() != "completed") {
        throw std::runtime_error(pj.string() + ": the run is " +
                                 (status != nullptr ? status->string() : std::string("unrecorded")) +
                                 "; --from needs a completed run");
    }
    const auto* refinement = process.find("refinement");
    if (refinement == nullptr || refinement->is_null()) {
        throw std::runtime_error(pj.string() + ": no refined geometry (PROCESS did not run)");
    }
    const fs::path ij = prev.directory / "input.json";
    const Json input = read_json_file(ij);
    if (const auto* x = process.find("experiment"); x != nullptr && !x->is_null()) {
        prev.model = read_final_model(*x, input, pj);
    }
    if (const auto* tt = refinement->find("detector_tilt_twist_mrad");
        tt != nullptr && tt->find("tilt") != nullptr && tt->find("twist") != nullptr) {
        prev.tilt_twist_mrad = std::pair{tt->find("tilt")->number(), tt->find("twist")->number()};
    } else if (const auto* v1 = refinement->find("tilt_twist_mrad"); v1 != nullptr && v1->items().size() == 2) {
        prev.tilt_twist_mrad = std::pair{v1->items()[0].number(), v1->items()[1].number()};  // version 1
    }
    for (const auto& o : process.find("outputs") != nullptr ? process.find("outputs")->items() : Json::Array{}) {
        const auto* kind = o.find("kind");
        const auto* p = o.find("path");
        if (kind != nullptr && p != nullptr && kind->string() == "gxparm") {
            prev.gxparm = fs::path(p->string()).is_absolute() ? fs::path(p->string())
                                                               : prev.directory / p->string();
        }
    }
    if (prev.gxparm.empty()) {
        prev.gxparm = prev.directory / "GXPARM.XDS";
    }
    // Without a final model (version 1) the orientation comes from GXPARM.XDS.
    if (!prev.model && !fs::exists(prev.gxparm)) {
        throw std::runtime_error("refined geometry not found: " + prev.gxparm.string());
    }

    const auto* inputs = input.find("inputs");
    const auto* images = inputs != nullptr ? inputs->find("images") : nullptr;
    if (images == nullptr || images->items().empty() || images->items()[0].find("path") == nullptr) {
        throw std::runtime_error(ij.string() + ": no image recorded");
    }
    prev.image = images->items()[0].find("path")->string();
    const auto* experiment = input.find("experiment");
    const auto* range = experiment != nullptr ? experiment->find("frame_range") : nullptr;
    const auto* value = range != nullptr ? range->find("value") : nullptr;
    if (value == nullptr || value->find("first") == nullptr || value->find("last") == nullptr) {
        throw std::runtime_error(ij.string() + ": no frame range recorded");
    }
    prev.frame_range = {static_cast<int>(value->find("first")->number()),
                        static_cast<int>(value->find("last")->number())};
    return prev;
}

}  // namespace mosflm::io
