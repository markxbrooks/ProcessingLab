/**
 * ExperimentModel: per-field metadata precedence and provenance for
 * ``ProcessingLab run`` (image header < XDS.INP < .xinfo < CLI), conflict
 * warnings, .xinfo sweep selection, and fresh output directories.
 */
#include "mosflm/io/ExperimentModel.hpp"
#include "mosflm/io/XinfoImporter.hpp"
#include "test_util.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

using mosflm::io::ExperimentInputs;
using mosflm::io::MetadataSource;
using mosflm::io::XinfoImporter;
using mosflm::io::create_output_directory;
using mosflm::io::resolve_experiment;
namespace fs = std::filesystem;

namespace {

void write(const fs::path& p, const std::string& text) {
    std::ofstream(p) << text;
}

bool has_warning(const mosflm::io::ExperimentModel& m, const std::string& what) {
    return std::any_of(m.warnings.begin(), m.warnings.end(), [&](const auto& w) {
        return w.rfind(what, 0) == 0;
    });
}

template <class F>
bool throws(F&& f) {
    try {
        f();
    } catch (const std::runtime_error&) {
        return true;
    }
    return false;
}

}  // namespace

int main() {
    const fs::path cbf = "testdata/cbf/sample_minicbf.cbf";
    CHECK_MSG(fs::exists(cbf), "missing sample_minicbf.cbf");
    const fs::path dir = fs::temp_directory_path() / "pl_experiment_model_test";
    fs::remove_all(dir);
    fs::create_directories(dir);
    for (const char* n : {"x_0001.cbf", "x_0002.cbf", "x_0003.cbf"}) {
        fs::copy_file(cbf, dir / n);
    }

    // Image only: everything from the header; the frames on disk set the range.
    ExperimentInputs in;
    in.image = dir / "x_0001.cbf";
    const auto a = resolve_experiment(in);
    CHECK(!a.hdf5 && a.image_template == "x_####.cbf");
    CHECK(a.wavelength_A.source == MetadataSource::ImageHeader);
    CHECK_CLOSE(a.wavelength_A.value, 0.9790, 1e-6);
    CHECK_CLOSE(a.distance_mm.value, 160.0, 1e-6);
    CHECK_EQ(a.frames_available, 3);
    CHECK(a.frame_range.value == std::pair(1, 3));
    CHECK(!a.cell.set() && a.warnings.empty());
    const auto lines = a.control_lines();
    CHECK(lines.size() == 2u && lines[1] == "TEMPLATE x_####.cbf");

    // XDS.INP overrides per field and warns where it disagrees with the header.
    write(dir / "XDS.INP",
          "DETECTOR_DISTANCE=165.0\nX-RAY_WAVELENGTH=0.979\n"
          "UNIT_CELL_CONSTANTS= 61.7 185.4 71.9 90 100.2 90\nSPACE_GROUP_NUMBER=3\n");
    in.xds_inp = dir / "XDS.INP";
    const auto b = resolve_experiment(in);
    CHECK(b.distance_mm.source == MetadataSource::XdsInp);
    CHECK_CLOSE(b.distance_mm.value, 165.0, 1e-6);
    CHECK(has_warning(b, "Distance: XDS.INP 165.00 mm overrides image header 160.00 mm"));
    CHECK(!has_warning(b, "Wavelength"));  // same value: no conflict
    CHECK(b.cell.source == MetadataSource::XdsInp && b.symmetry.value == 3);
    CHECK(b.pixel_mm.source == MetadataSource::ImageHeader);  // XDS.INP has none

    // .xinfo outranks XDS.INP, but only for the fields it has: its lack of
    // detector geometry must not erase the distance.
    write(dir / "one.xinfo",
          "BEGIN WAVELENGTH NATIVE\nWAVELENGTH 0.9800\nEND WAVELENGTH NATIVE\n"
          "BEGIN SWEEP S1\nWAVELENGTH NATIVE\nSTART_END 2 3\nEND SWEEP S1\n");
    in.xinfo = dir / "one.xinfo";
    const auto c = resolve_experiment(in);
    CHECK(c.wavelength_A.source == MetadataSource::Xinfo);
    CHECK(has_warning(c, "Wavelength: .xinfo 0.98000"));
    CHECK(c.distance_mm.source == MetadataSource::XdsInp);
    CHECK(c.frame_range.source == MetadataSource::Xinfo && c.frame_range.value == std::pair(2, 3));

    // The command line outranks everything; out-of-range frames are refused.
    in.range = std::pair(1, 2);
    CHECK(resolve_experiment(in).frame_range.source == MetadataSource::CommandLine);
    in.range = std::pair(1, 9);
    CHECK(throws([&] { (void)resolve_experiment(in); }));
    in.range.reset();

    // Several sweeps: refuse to guess; --sweep picks one and its wavelength.
    write(dir / "two.xinfo",
          "BEGIN WAVELENGTH W1\nWAVELENGTH 0.9790\nEND WAVELENGTH W1\n"
          "BEGIN WAVELENGTH W2\nWAVELENGTH 1.2000\nEND WAVELENGTH W2\n"
          "BEGIN SWEEP A\nWAVELENGTH W1\nSTART_END 1 2\nEND SWEEP A\n"
          "BEGIN SWEEP B\nWAVELENGTH W2\nSTART_END 2 3\nEND SWEEP B\n");
    CHECK(XinfoImporter::sweep_names_in_file(dir / "two.xinfo") == std::vector<std::string>({"A", "B"}));
    in.xinfo = dir / "two.xinfo";
    CHECK(throws([&] { (void)resolve_experiment(in); }));
    in.sweep = "B";
    const auto e = resolve_experiment(in);
    CHECK_CLOSE(e.wavelength_A.value, 1.2, 1e-9);
    CHECK(e.frame_range.value == std::pair(2, 3));
    in.sweep = "Z";
    CHECK(throws([&] { (void)resolve_experiment(in); }));

    // No image at all.
    CHECK(throws([] { (void)resolve_experiment(ExperimentInputs{}); }));

    // Output directories: an explicit one must be new; default names never
    // collide, even within the same second.
    const auto out = create_output_directory(dir / "results" / "run1", dir, 0);
    CHECK(fs::is_directory(out));
    CHECK(throws([&] { (void)create_output_directory(dir / "results" / "run1", dir, 0); }));
    const std::time_t t = 1792700000;
    const auto d1 = create_output_directory(std::nullopt, dir, t);
    const auto d2 = create_output_directory(std::nullopt, dir, t);
    CHECK(d1 != d2 && fs::is_directory(d1) && fs::is_directory(d2));
    CHECK(d1.filename().string().rfind("processing_", 0) == 0);
    CHECK(d2.filename().string() == d1.filename().string() + "_2");

    std::error_code ec;
    fs::remove_all(dir, ec);
    return mosflm_test::finish("experiment_model");
}
