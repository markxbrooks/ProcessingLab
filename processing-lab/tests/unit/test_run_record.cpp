/**
 * Run records: Json serialisation, atomic writes, planned stages from
 * CONTROL, and the process.json lifecycle (running → completed / failed /
 * skipped) with process.md rendered from the records.
 */
#include "mosflm/io/Json.hpp"
#include "mosflm/io/RunRecord.hpp"
#include "test_util.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
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

    return mosflm_test::finish("run_record");
}
