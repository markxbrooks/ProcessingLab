/**
 * Run records: Json serialisation, atomic writes, planned stages from
 * CONTROL, and the process.json lifecycle (running → completed / failed /
 * skipped) with process.md rendered from the records.
 */
#include "mosflm/geometry/LabFrame.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/io/Json.hpp"
#include "mosflm/io/RunRecord.hpp"
#include "test_util.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using mosflm::io::Json;
using mosflm::io::ProcessRecord;
namespace fs = std::filesystem;

namespace {

std::string slurp(const fs::path& p) {
    std::ifstream in(p);
    std::ostringstream s;
    s << in.rdbuf();
    return s.str();
}

bool contains(const std::string& text, const std::string& what) {
    return text.find(what) != std::string::npos;
}

}  // namespace

int main() {
    // Json: ordered objects, escaping, integers vs reals, non-finite → null.
    Json j = Json::object();
    j["text"] = "a \"quoted\"\nline\t\x01";
    j["int"] = 3;
    j["real"] = 0.25;
    j["nan"] = std::numeric_limits<double>::quiet_NaN();
    j["cell"] = Json(Json::Array{Json(78.19), Json(90)});
    j["empty"] = Json::array();
    const std::string d = j.dump();
    CHECK(contains(d, "\"text\": \"a \\\"quoted\\\"\\nline\\t\\u0001\""));
    CHECK(contains(d, "\"int\": 3,"));
    CHECK(contains(d, "\"real\": 0.25,"));
    CHECK(contains(d, "\"nan\": null,"));
    CHECK(contains(d, "\"cell\": [78.19, 90],"));  // numbers stay on one line
    CHECK(d.find("\"text\"") < d.find("\"int\""));  // insertion order kept
    CHECK(j.find("int") != nullptr && j.find("int")->number() == 3.0);
    CHECK(j.find("missing") == nullptr);

    // Parse: what dump writes reads back the same; escapes and errors.
    CHECK(Json::parse(d).dump() == j.dump());
    {
        const Json p = Json::parse(R"({"s": "caf\u00e9 \"q\" \ud83d\ude00", "n": -1.5e2, "i": 42,
                                      "a": [true, false, null], "o": {}})");
        CHECK(p.find("s")->string() == "caf\xc3\xa9 \"q\" \xf0\x9f\x98\x80");
        CHECK(p.find("n")->number() == -150.0);
        CHECK(p.find("i")->number() == 42.0);
        CHECK(p.find("a")->items().size() == 3 && p.find("a")->items()[0].boolean());
        CHECK(p.find("a")->items()[2].is_null() && p.find("o")->is_object());
        const auto bad = [](const char* text) {
            try {
                (void)Json::parse(text);
            } catch (const std::runtime_error&) {
                return true;
            }
            return false;
        };
        CHECK(bad("{\"a\": 1,}"));
        CHECK(bad("[1 2]"));
        CHECK(bad("{\"a\": 1} x"));
        CHECK(bad("\"open"));
    }

    // Atomic write: content in place, no temporary left.
    const fs::path dir = mosflm_test::output_dir("run_record");
    mosflm::io::write_file_atomic(dir / "x.json", j.dump());
    CHECK(slurp(dir / "x.json") == d);
    CHECK(!fs::exists(dir / "x.json.tmp"));

    // Planned stages: GO runs the first pending step (FINDSPOTS, AUTOINDEX,
    // PROCESS, STRATEGY); nothing after EXIT.
    using mosflm::io::planned_stages;
    CHECK(planned_stages("FINDSPOTS FIND 1\nGO\nAUTOINDEX DPS\nGO\nPROCESS 1 TO 5\nGO\nEXIT\nSTRATEGY\nGO\n") ==
          std::vector<std::string>({"findspots", "autoindex", "process"}));
    CHECK(planned_stages("AUTOINDEX\nFINDSPOTS\n! comment\nGO\nrun\n") ==
          std::vector<std::string>({"findspots", "autoindex"}));
    CHECK(planned_stages("DISTANCE 100\nEXIT\n").empty());

    // Lifecycle: created running, stages pending; findspots completes,
    // autoindex fails → run failed with the error, in process.md too.
    Json input = Json::object();
    input["mode"] = "run";
    input["warnings"] = Json(Json::Array{Json("Distance: XDS.INP 165 overrides image header 160")});
    {
        ProcessRecord rec(dir, input, {"findspots", "autoindex", "process"});
        CHECK(fs::exists(dir / "input.json") && fs::exists(dir / "process.md"));
        std::string pj = slurp(dir / "process.json");
        CHECK(contains(pj, "\"status\": \"running\""));
        CHECK(contains(pj, "\"status\": \"pending\""));
        CHECK(contains(pj, "Distance: XDS.INP 165"));  // input warnings carried over
        // Version 2: every stage key present from the start; validation is
        // separate from completion; the run directory is recorded.
        CHECK(contains(pj, "\"schema_version\": 2"));
        CHECK(contains(pj, "\"started_at\": \"20"));
        CHECK(contains(pj, "\"validation\": {\n    \"status\": \"not_assessed\""));
        CHECK(contains(pj, "\"run_directory\": \"" + fs::absolute(dir).lexically_normal().string() + "\""));
        {
            const Json doc = Json::parse(pj);
            const Json& first = doc.find("stages")->items()[0];
            for (const char* key : {"name", "status", "started_at", "finished_at", "duration_s", "message",
                                    "results", "error"}) {
                CHECK_MSG(first.find(key) != nullptr, std::string("stage key ") + key);
            }
            // RFC 3339: the offset has a colon (…-07:00 or …+00:00).
            const std::string at = doc.find("started_at")->string();
            CHECK_MSG(at.size() == 25 && at[22] == ':', "RFC 3339 offset: " + at);
        }

        mosflm::control::ControlState st;
        st.findspots_pending = true;
        rec.begin_stage(st);
        CHECK(contains(slurp(dir / "process.json"), "\"name\": \"findspots\",\n      \"status\": \"running\""));
        mosflm::control::ControlGoResult r;
        r.ran_findspots = true;
        st.findspots_pending = false;
        st.findspots_n_spots = 42;
        rec.end_stage(r, st);
        pj = slurp(dir / "process.json");
        CHECK(contains(pj, "\"name\": \"findspots\",\n      \"status\": \"completed\""));
        CHECK(contains(pj, "\"spots\": 42"));
        CHECK(Json::parse(pj).find("stages")->items()[0].find("duration_s")->is_number());

        st.autoindex_pending = true;
        rec.begin_stage(st);
        rec.fail("indexing failed: too few spots");
        pj = slurp(dir / "process.json");
        CHECK(contains(pj, "\"name\": \"autoindex\",\n      \"status\": \"failed\""));
        CHECK(contains(pj, "\"error\": \"indexing failed: too few spots\""));
        CHECK(contains(pj, "\"name\": \"process\",\n      \"status\": \"pending\""));
        const std::string md = slurp(dir / "process.md");
        CHECK(contains(md, "**Status:** failed"));
        CHECK(contains(md, "> **Error:** indexing failed: too few spots"));
        CHECK(contains(md, "| findspots | completed | 42 spots |"));
    }

    // Completing with stages never reached marks them skipped.
    {
        const fs::path dir2 = mosflm_test::output_dir("run_record_skip");
        ProcessRecord rec(dir2, input, {"findspots", "process"});
        rec.complete();
        const std::string pj = slurp(dir2 / "process.json");
        CHECK(contains(pj, "\"status\": \"completed\""));
        CHECK(contains(pj, "\"status\": \"skipped\""));
    }

    // run --from: a completed run's image, range and refined GXPARM.XDS.
    {
        using mosflm::io::read_previous_run;
        const fs::path prev = mosflm_test::output_dir("run_record_previous");
        const auto refused = [](const fs::path& p, const std::string& what) {
            try {
                (void)read_previous_run(p);
            } catch (const std::runtime_error& ex) {
                return std::string(ex.what()).find(what) != std::string::npos;
            }
            return false;
        };
        CHECK(refused(prev, "no process.json"));
        mosflm::io::write_file_atomic(prev / "input.json", R"({"inputs": {"images": [{"path": "/data/x_0001.cbf"}]},
  "experiment": {"frame_range": {"value": {"first": 2, "last": 90}}}})");
        mosflm::io::write_file_atomic(prev / "process.json", R"({"status": "failed", "refinement": null})");
        CHECK(refused(prev, "the run is failed"));
        mosflm::io::write_file_atomic(prev / "process.json", R"({"status": "completed", "refinement": null})");
        CHECK(refused(prev, "no refined geometry"));
        mosflm::io::write_file_atomic(prev / "process.json",
                                      R"({"status": "completed", "refinement": {"tilt_twist_mrad": [9.4, -8.4]},
  "outputs": [{"kind": "gxparm", "path": "GXPARM.XDS"}]})");
        CHECK(refused(prev, "refined geometry not found"));
        mosflm::io::write_file_atomic(prev / "GXPARM.XDS", "! placeholder\n");
        const auto pr = read_previous_run(prev / "process.json");  // file or directory
        CHECK(pr.directory == fs::absolute(prev).lexically_normal());
        CHECK(pr.image == fs::path("/data/x_0001.cbf"));
        CHECK(pr.frame_range == std::pair(2, 90));
        CHECK(pr.gxparm == pr.directory / "GXPARM.XDS");
        CHECK(pr.tilt_twist_mrad == std::pair(9.4, -8.4));  // version 1 layout

        // Version 2 layout; a newer schema is refused.
        mosflm::io::write_file_atomic(prev / "process.json",
                                      R"({"schema_version": 2, "status": "completed",
  "refinement": {"detector_tilt_twist_mrad": {"tilt": 1.5, "twist": -2.5}},
  "outputs": [{"kind": "gxparm", "path": "GXPARM.XDS"}]})");
        CHECK(read_previous_run(prev).tilt_twist_mrad == std::pair(1.5, -2.5));
        mosflm::io::write_file_atomic(prev / "process.json", R"({"schema_version": 3, "status": "completed"})");
        CHECK(refused(prev, "schema_version 3 is newer"));

        // Schema 2 final model: carried forward only if it is complete and
        // self-consistent (matrix ↔ cell, cell ↔ space group, λ / pixel ↔
        // input.json).
        const mosflm::UnitCell tet{57.938, 57.938, 150.244, 90.0, 90.0, 90.0};
        const double lam = 0.9795;
        const Eigen::Matrix3d rot =
            Eigen::AngleAxisd(0.7, Eigen::Vector3d(0.3, -0.5, 0.8).normalized()).toRotationMatrix();
        const auto amat_for = [&](const mosflm::UnitCell& c) {
            return mosflm::cpp_ub_to_mosflm_amat(mosflm::OrientationMatrix::from_cell_and_rotation(c, rot).ub(), lam);
        };
        const auto model_json = [&](const Eigen::Matrix3d& a, double wavelength, const char* convention,
                                    int sg, double b_axis) {
            std::ostringstream s;
            s << std::setprecision(12) << R"({"schema_version": 2, "status": "completed",
  "refinement": {"stage": "process"},
  "experiment": {"crystal": {"unit_cell": {"a": 57.938, "b": )" << b_axis << R"(, "c": 150.244,
      "alpha": 90, "beta": 90, "gamma": 90}, "space_group": {"number": )" << sg << R"(}},
    "wavelength_angstrom": )" << wavelength << R"(,
    "detector": {"distance_mm": 199.98, "beam_centre_mm": {"x": 166.12, "y": 172.93},
      "pixel_size_mm": {"x": 0.075, "y": 0.075}, "tilt_twist_mrad": {"tilt": 9.39, "twist": -8.43}},
    "mosaicity_deg": 0.101,
    "orientation": {"convention": ")" << convention << R"(", "wavelength_angstrom": )" << lam
              << R"(, "invertx": false, "matrix": [)";
            for (int r = 0; r < 3; ++r) {
                s << (r ? ", " : "") << '[' << a(r, 0) << ", " << a(r, 1) << ", " << a(r, 2) << ']';
            }
            s << "]}}}";
            return s.str();
        };
        mosflm::io::write_file_atomic(prev / "input.json", R"({"inputs": {"images": [{"path": "/data/t_master.h5"}]},
  "experiment": {"frame_range": {"value": {"first": 1, "last": 100}},
    "wavelength": {"value": 0.9795}, "pixel_size": {"value": [0.075, 0.075]}}})");
        fs::remove(prev / "GXPARM.XDS");  // not needed with a final model
        mosflm::io::write_file_atomic(prev / "process.json", model_json(amat_for(tet), lam, "mosflm_amat", 75, 57.938));
        const auto good = read_previous_run(prev);
        CHECK(good.model.has_value());
        CHECK(good.model->space_group == 75);
        CHECK(std::fabs(good.model->cell[2] - 150.244) < 1e-9);
        CHECK(good.model->tilt_twist_mrad == std::pair(9.39, -8.43));
        CHECK(std::fabs(good.model->mosaicity_deg - 0.101) < 1e-12);
        CHECK(!good.model->invertx);
        CHECK(std::fabs(good.model->matrix_wavelength_A - lam) < 1e-12);
        // A matrix built for another cell (the indexing estimate) is refused.
        const mosflm::UnitCell indexing{57.94, 58.24, 150.33, 90.45, 90.18, 90.02};
        mosflm::io::write_file_atomic(prev / "process.json",
                                      model_json(amat_for(indexing), lam, "mosflm_amat", 75, 57.938));
        CHECK(refused(prev, "the orientation matrix implies cell"));
        // λ disagreeing with input.json, a cell the space group forbids, an
        // unknown convention, and missing geometry are refused too.
        mosflm::io::write_file_atomic(prev / "process.json", model_json(amat_for(tet), 1.0, "mosflm_amat", 75, 57.938));
        CHECK(refused(prev, "final wavelength 1.000000 differs from 0.979500"));
        mosflm::io::write_file_atomic(prev / "process.json", model_json(amat_for(tet), lam, "mosflm_amat", 75, 58.5));
        CHECK(refused(prev, "does not fit space group 75"));
        mosflm::io::write_file_atomic(prev / "process.json", model_json(amat_for(tet), lam, "xds", 75, 57.938));
        CHECK(refused(prev, "orientation convention 'xds' is not mosflm_amat"));
        mosflm::io::write_file_atomic(prev / "process.json", R"({"schema_version": 2, "status": "completed",
  "refinement": {}, "experiment": {"crystal": {}}})");
        CHECK(refused(prev, "experiment.unit_cell is missing"));

        // A multi-sweep run names its sweep directories.
        const fs::path multi = mosflm_test::output_dir("run_record_multi");
        fs::create_directories(multi / "sweep_1");
        mosflm::io::write_file_atomic(multi / "sweep_1" / "process.json", "{}");
        CHECK(refused(multi, "holds several sweeps (sweep_1)"));
    }

    // --run-mode: names, parsing, and inputs that must match the mode.
    {
        using mosflm::io::RunMode;
        using mosflm::io::parse_run_mode;
        using mosflm::io::validate_run_mode;
        CHECK(parse_run_mode("xds") == RunMode::XDS);
        CHECK(parse_run_mode("XIA2") == RunMode::XIA2);
        CHECK(parse_run_mode("MosFlm") == RunMode::MOSFLM);
        CHECK(parse_run_mode("processinglab") == RunMode::PROCESSINGLAB);
        CHECK(!parse_run_mode("dials").has_value());
        CHECK(std::string(mosflm::io::run_mode_name(RunMode::XDS)) == "RunMode.XDS");
        CHECK(std::string(mosflm::io::run_mode_name(RunMode::MOSFLM)) == "RunMode.MOSFLM");
        CHECK(std::string(mosflm::io::run_mode_name(RunMode::XIA2)) == "RunMode.XIA2");
        const auto refused = [](RunMode m, bool script, bool xds, bool xinfo, bool sweep) {
            try {
                validate_run_mode(m, script, xds, xinfo, sweep);
            } catch (const std::invalid_argument&) {
                return true;
            }
            return false;
        };
        CHECK(!refused(RunMode::PROCESSINGLAB, false, true, true, true));  // any imports
        CHECK(!refused(RunMode::XDS, false, true, false, false));
        CHECK(refused(RunMode::XDS, false, false, false, false));   // needs --xds-inp
        CHECK(refused(RunMode::XDS, false, true, true, false));     // not with --xinfo
        CHECK(!refused(RunMode::XIA2, false, false, true, true));
        CHECK(refused(RunMode::XIA2, false, false, false, false));  // needs --xinfo
        CHECK(refused(RunMode::XIA2, false, true, true, false));    // not with --xds-inp
        CHECK(!refused(RunMode::MOSFLM, true, false, false, false));
        CHECK(refused(RunMode::MOSFLM, false, false, false, false));  // not with 'run'
        CHECK(refused(RunMode::XDS, true, true, false, false));       // not a script run
        // --xparm / --from: an import like --xds-inp, refused where the mode
        // names a single source.
        const auto refused_xparm = [](RunMode m) {
            try {
                validate_run_mode(m, false, m == RunMode::XDS, m == RunMode::XIA2, false, true);
            } catch (const std::invalid_argument&) {
                return true;
            }
            return false;
        };
        CHECK(!refused_xparm(RunMode::PROCESSINGLAB));
        CHECK(refused_xparm(RunMode::XDS));
        CHECK(refused_xparm(RunMode::XIA2));
    }

    return mosflm_test::finish("run_record");
}
