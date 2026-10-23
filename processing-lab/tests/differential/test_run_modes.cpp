/**
 * End-to-end input modes of the ``ProcessingLab`` binary on the CCP4 insulin
 * frame, with fixtures from ``testdata/run_modes``:
 *
 *   - ``run --image`` (header metadata, default timestamped output directory)
 *   - ``run --image --xds-inp`` (per-field provenance from XDS.INP)
 *   - ``run --xinfo`` (image named by the import)
 *   - ``run --xinfo`` with two sweeps (refused without ``--sweep``)
 *   - ``run --image --image`` (one ``sweep_N`` directory per image)
 *   - ``run`` refusals: existing output directory, frames not on disk
 *   - CONTROL script mode with ``--output-dir``
 *   - run records: input.json provenance, process.json lifecycle (completed
 *     and failed), process.md
 *
 * Every run happens in ``<build>/test_output/run_modes``; the test checks that
 * nothing is written into ``testdata/run_modes`` or next to the image.
 * Soft-skips if the insulin image is missing.
 */
#include "test_util.hpp"

#include <gemmi/mtz.hpp>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <sys/wait.h>

namespace fs = std::filesystem;

namespace {

fs::path insulin_image() {
    if (const char* env = std::getenv("MOSFLM_IMAGE")) {
        return env;
    }
    const char* ccp4 = std::getenv("CCP4");
    return fs::path(ccp4 && *ccp4 ? ccp4 : "/Applications/ccp4-9") /
           "examples/data/insulin_1_001.img";
}

std::string slurp(const fs::path& p) {
    std::ifstream in(p);
    std::ostringstream s;
    s << in.rdbuf();
    return s.str();
}

/// Copy ``src`` to ``dst`` with ``@IMAGE_DIR@`` filled in.
void render(const fs::path& src, const fs::path& dst, const fs::path& image_dir) {
    std::string text = slurp(src);
    const std::string key = "@IMAGE_DIR@";
    for (auto pos = text.find(key); pos != std::string::npos; pos = text.find(key, pos)) {
        text.replace(pos, key.size(), image_dir.string());
    }
    std::ofstream(dst) << text;
}

void write_text_file(const fs::path& p, const std::string& text) {
    std::ofstream(p) << text;
}

/// Run ``ProcessingLab args`` in ``cwd``; stdout+stderr to ``cwd/<log>``.
int run(const fs::path& cwd, const std::string& args, const std::string& log) {
    const std::string cmd = "cd '" + cwd.string() + "' && '" + std::string(PROCESSINGLAB_BIN) +
                            "' " + args + " > " + log + " 2>&1";
    const int rc = std::system(cmd.c_str());
    return WIFEXITED(rc) ? WEXITSTATUS(rc) : -1;
}

/// Directory listing with sizes and mtimes, to detect stray writes.
std::map<std::string, std::pair<std::uintmax_t, fs::file_time_type>> snapshot(const fs::path& dir) {
    std::map<std::string, std::pair<std::uintmax_t, fs::file_time_type>> out;
    for (const auto& e : fs::directory_iterator(dir)) {
        if (e.is_regular_file()) {
            out[e.path().filename().string()] = {e.file_size(), e.last_write_time()};
        } else {
            out[e.path().filename().string() + "/"] = {0, fs::file_time_type{}};
        }
    }
    return out;
}

bool contains(const std::string& text, const std::string& what) {
    return text.find(what) != std::string::npos;
}

/// The six standard outputs, with an MTZ of usable size.
void check_outputs(const fs::path& dir, const std::string& spots, const std::string& label) {
    for (const std::string& f : std::vector<std::string>{"process.log", "process.mtz", "XPARM.XDS",
                                                         "GXPARM.XDS", "ProcessingLabIndex.txt", spots,
                                                         "input.json", "process.json", "process.md"}) {
        CHECK_MSG(fs::exists(dir / f), label + ": missing " + f);
    }
    // Lifecycle finished, no temporary files left by the atomic writes.
    CHECK_MSG(contains(slurp(dir / "process.json"), "\"status\": \"completed\""),
              label + ": process.json not completed");
    CHECK_MSG(!fs::exists(dir / "process.json.tmp"), label + ": stray process.json.tmp");
    if (fs::exists(dir / "process.mtz")) {
        gemmi::Mtz mtz;
        mtz.read_file((dir / "process.mtz").string());
        CHECK_MSG(mtz.nreflections > 500, label + ": MTZ has too few reflections");
    }
}

}  // namespace

int main() {
    const fs::path image = insulin_image();
    if (!fs::exists(image)) {
        std::cout << "run_modes: SKIP (no insulin image; set MOSFLM_IMAGE or CCP4)\n";
        return 0;
    }
    const fs::path fixtures = fs::absolute("testdata/run_modes");
    CHECK_MSG(fs::exists(fixtures / "process.in"), "missing testdata/run_modes");
    const fs::path work = mosflm_test::output_dir("run_modes");
    const fs::path image_dir = fs::absolute(image).parent_path();
    for (const char* f : {"insulin_XDS.INP", "insulin.xinfo", "insulin_two_sweeps.xinfo"}) {
        render(fixtures / f, work / f, image_dir);
    }
    const auto fixtures_before = snapshot(fixtures);
    const auto image_dir_before = snapshot(image_dir);
    const std::string img = "'" + image.string() + "'";

    // 1. run --image: header metadata, default processing_YYYYMMDD_HHMMSS.
    {
        CHECK_EQ(run(work, "run --image " + img, "image_only.out"), 0);
        // The run reports where it wrote (counting processing_* directories
        // would race another test run in the same build tree).
        const std::string out = slurp(work / "image_only.out");
        const std::string key = "Output directory: ";
        const auto pos = out.rfind(key);
        CHECK_MSG(pos != std::string::npos, "run did not report its output directory");
        if (pos != std::string::npos) {
            const fs::path made = out.substr(pos + key.size(), out.find('\n', pos) - pos - key.size());
            CHECK(made.parent_path() == work);
            CHECK(made.filename().string().rfind("processing_", 0) == 0);
            check_outputs(made, "spots.spt", "run --image");
            const std::string log = slurp(made / "process.log");
            CHECK(contains(log, "Wavelength:    0.97900 Å  [image header]"));
            CHECK(contains(log, "Cell:          from autoindexing"));
        }
    }

    // 2. run --image --xds-inp: cell and symmetry from XDS.INP.
    {
        CHECK_EQ(run(work, "run --image " + img + " --xds-inp insulin_XDS.INP --output-dir xds",
                     "xds.out"), 0);
        check_outputs(work / "xds", "spots.spt", "run --xds-inp");
        const std::string log = slurp(work / "xds" / "process.log");
        CHECK(contains(log, "Cell:          78.19 78.19 78.19 90.00 90.00 90.00  [XDS.INP]"));
        CHECK(contains(log, "Space group:   197  [XDS.INP]"));
        CHECK(contains(log, "Image source:  " + image.string() + "  [CLI]"));
    }

    // 3. run --xinfo with no --image: the import names the image.
    {
        CHECK_EQ(run(work, "run --xinfo insulin.xinfo --output-dir xinfo", "xinfo.out"), 0);
        check_outputs(work / "xinfo", "spots.spt", "run --xinfo");
        const std::string log = slurp(work / "xinfo" / "process.log");
        CHECK(contains(log, "[.xinfo]"));
        CHECK(contains(log, "Image source:  " + (image_dir / "insulin_1_001.img").string()));
    }

    // 4. Two sweeps: refused without --sweep (and nothing created); --sweep
    //    picks one, with the wavelength that sweep names.
    {
        CHECK_EQ(run(work, "run --xinfo insulin_two_sweeps.xinfo --output-dir ambiguous",
                     "ambiguous.out"), 2);
        CHECK(contains(slurp(work / "ambiguous.out"), "choose one with --sweep"));
        CHECK(!fs::exists(work / "ambiguous"));
        CHECK_EQ(run(work, "run --xinfo insulin_two_sweeps.xinfo --sweep SWEEP2 --output-dir sweep2",
                     "sweep2.out"), 0);
        CHECK(contains(slurp(work / "sweep2" / "process.log"), "Wavelength:    0.98500 Å  [.xinfo]"));
    }

    // 5. Two images: one sweep_N directory each.
    {
        CHECK_EQ(run(work, "run --image " + img + " --image " + img + " --output-dir multi",
                     "multi.out"), 0);
        check_outputs(work / "multi" / "sweep_1", "spots.spt", "sweep_1");
        check_outputs(work / "multi" / "sweep_2", "spots.spt", "sweep_2");
    }

    // 6. Refusals: an existing output directory is never reused; frames must
    //    exist. Exit code 2, a clear message, nothing created.
    {
        CHECK_EQ(run(work, "run --image " + img + " --output-dir xds", "reuse.out"), 2);
        CHECK(contains(slurp(work / "reuse.out"), "already exists"));
        CHECK_EQ(run(work, "run --image " + img + " --range 1:5 --output-dir badrange",
                     "badrange.out"), 2);
        CHECK(contains(slurp(work / "badrange.out"), "outside frames 1–1"));
        CHECK(!fs::exists(work / "badrange"));
    }

    // 7. CONTROL script mode with --output-dir.
    {
        CHECK_EQ(run(work, "--image " + img + " --output-dir script '" +
                     (fixtures / "process.in").string() + "'", "script.out"), 0);
        check_outputs(work / "script", "insulin.spt", "script mode");
        CHECK(contains(slurp(work / "script" / "process.log"), "***** CONTROL INPUT *****"));
    }

    // 8. A failing stage: process.json ends "failed" with the error, on the
    //    stage that failed; later stages stay pending.
    {
        write_text_file(work / "fail.in", "AUTOINDEX DPS THRESHOLD 15\nGO\nPROCESS 1 TO 1\nGO\nEXIT\n");
        CHECK_EQ(run(work, "--image " + img + " --output-dir failed fail.in", "fail.out"), 1);
        const std::string pj = slurp(work / "failed" / "process.json");
        CHECK(contains(pj, "\"status\": \"failed\""));
        CHECK(contains(pj, "\"error\": \"SpotList: cannot open"));
        CHECK(contains(pj, "\"name\": \"process\",\n      \"status\": \"pending\""));
        CHECK(contains(slurp(work / "failed" / "process.md"), "**Status:** failed"));
    }

    // 9. input.json provenance: CONTROL overrides the header, history kept.
    {
        const std::string ij = slurp(work / "xds" / "input.json");
        CHECK(contains(ij, "\"source\": \"xds_inp\""));
        CHECK(contains(ij, "\"overrides\": [\n        {\n          \"value\": 160,\n          \"source\": \"image_header\""));
    }

    // Nothing written into the fixtures or next to the image.
    CHECK_MSG(snapshot(fixtures) == fixtures_before, "testdata/run_modes was modified");
    CHECK_MSG(snapshot(image_dir) == image_dir_before, "files written next to the image");

    return mosflm_test::finish("run_modes");
}
