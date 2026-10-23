/**
 * Unit tests for Xia2 .xinfo and XDS.INP importers.
 */

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/io/ExperimentModel.hpp"
#include "mosflm/io/ImportedSweepGeometry.hpp"
#include "mosflm/io/XdsInpImporter.hpp"
#include "mosflm/io/XinfoImporter.hpp"
#include "test_util.hpp"

#include <cmath>
#include <filesystem>
#include <string>

using mosflm::control::ControlSession;
using mosflm::io::ImportedSweepGeometry;
using mosflm::io::ExperimentModel;
using mosflm::io::MetadataSource;
using mosflm::io::XdsInpImporter;
using mosflm::io::XinfoImporter;
using mosflm::io::apply_imported_geometry;
using mosflm::io::write_xds_inp;
using mosflm::io::template_from_image_name;
using mosflm::io::xds_question_marks_to_hash;

int main() {
    {
        CHECK(
            xds_question_marks_to_hash("sweep_1_0???.cbf") == "sweep_1_0###.cbf"
        );
        CHECK(template_from_image_name("sweep_1_0001.cbf") == "sweep_1_####.cbf");
    }

    {
        const auto g = XinfoImporter::parse_file("testdata/import/sample.xinfo");
        CHECK(g.has_wavelength);
        CHECK(std::fabs(g.wavelength_A - 0.979490) < 1e-9);
        CHECK(g.has_directory);
        CHECK(g.directory.string() == "/data/sweep");
        CHECK(g.has_image);
        CHECK(g.image_path.string() == "sweep_1_0001.cbf");
        CHECK(g.has_process_range);
        CHECK_EQ(g.process_first, 1);
        CHECK_EQ(g.process_last, 900);
        CHECK(g.has_template);
        CHECK(g.image_template == "sweep_1_####.cbf");

        const auto lines = g.to_control_lines();
        bool saw_process = false;
        bool saw_template = false;
        bool saw_image = false;
        for (const auto& line : lines) {
            if (line.rfind("PROCESS ", 0) == 0) {
                saw_process = true;
                CHECK(line.find("1 TO 900") != std::string::npos);
            }
            if (line.rfind("TEMPLATE ", 0) == 0) {
                saw_template = true;
            }
            if (line.rfind("IMAGE ", 0) == 0) {
                saw_image = true;
            }
        }
        CHECK(saw_process);
        CHECK(saw_template);
        CHECK(!saw_image);  // multi-frame: TEMPLATE only
    }

    {
        const auto g =
            XdsInpImporter::parse_file("testdata/import/sample_xds.inp");
        CHECK(g.has_wavelength);
        CHECK(g.has_distance);
        CHECK(std::fabs(g.distance_mm - 385.190) < 1e-6);
        CHECK(g.has_pixel);
        CHECK(std::fabs(g.pixel_x_mm - 0.172) < 1e-9);
        CHECK(g.has_beam);
        CHECK(std::fabs(g.beam_x_mm - 1245.82 * 0.172) < 1e-4);
        CHECK(std::fabs(g.beam_y_mm - 1209.07 * 0.172) < 1e-4);
        CHECK(g.has_cell);
        CHECK(std::fabs(g.cell_a - 61.70) < 1e-6);
        CHECK(g.has_symmetry);
        CHECK_EQ(g.symmetry_number, 3);
        CHECK(g.has_process_range);
        CHECK_EQ(g.process_first, 1);
        CHECK_EQ(g.process_last, 900);
        CHECK(g.has_oscillation);
        CHECK(std::fabs(g.oscillation_deg - 0.20) < 1e-9);
        CHECK(g.has_template);
        CHECK(g.image_template == "sweep_1_0###.cbf");
        CHECK(g.has_directory);
        CHECK(g.directory.string() == "/data/results/old");

        bool saw_beam = false;
        bool saw_angle = false;
        for (const auto& line : g.to_control_lines()) {
            if (line.rfind("BEAM ", 0) == 0) {
                saw_beam = true;
            }
            if (line.find("ANGLE") != std::string::npos) {
                saw_angle = true;
            }
        }
        CHECK(saw_beam);
        CHECK(saw_angle);
    }

    {
        const auto g =
            XdsInpImporter::parse_file("testdata/import/sample_xds_no_org.inp");
        CHECK(!g.has_beam);
        CHECK(g.has_distance);
        CHECK(g.has_pixel);
        for (const auto& line : g.to_control_lines()) {
            CHECK_MSG(
                line.rfind("BEAM ", 0) != 0,
                "integrate-style INP must not force BEAM"
            );
        }
    }

    {
        // Merge order: XDS then xinfo — xinfo wins paths / range / λ.
        ControlSession session;
        session.set_output_directory(mosflm_test::output_dir("xinfo_xds_import_121"));
        apply_imported_geometry(
            session, XdsInpImporter::parse_file("testdata/import/sample_xds.inp")
        );
        apply_imported_geometry(
            session, XinfoImporter::parse_file("testdata/import/sample.xinfo")
        );
        const auto& st = session.state();
        CHECK(st.has_wavelength);
        CHECK(std::fabs(st.wavelength_A - 0.979490) < 1e-9);
        CHECK(st.has_distance);
        CHECK(st.has_beam);
        CHECK(st.has_cell);
        CHECK_EQ(st.symmetry_number, 3);
        CHECK(st.directory.string() == "/data/sweep");
        CHECK(st.image_template == "sweep_1_####.cbf");
        CHECK_EQ(st.process_first, 1);
        CHECK_EQ(st.process_last, 900);
        CHECK(st.process_pending);
        CHECK(st.has_process_phi);  // ANGLE from XDS PROCESS line
        CHECK(std::fabs(st.process_angle_deg - 0.20) < 1e-9);
    }

    {
        const auto g = XinfoImporter::parse_string(
            "BEGIN WAVELENGTH NATIVE\n"
            "WAVELENGTH 1.0\n"
            "END WAVELENGTH NATIVE\n"
            "BEGIN SWEEP S1\n"
            "DIRECTORY /tmp\n"
            "IMAGE one.cbf\n"
            "START_END 5 5\n"
            "END SWEEP S1\n"
        );
        CHECK(g.has_image);
        CHECK(!g.has_template);  // single frame
        bool saw_image = false;
        for (const auto& line : g.to_control_lines()) {
            if (line.rfind("IMAGE ", 0) == 0) {
                saw_image = true;
            }
        }
        CHECK(saw_image);
    }

    {
        // Resolved options written at the start of a run, readable again as
        // XDS.INP. ``#`` in the template is stored as ``?``.
        ExperimentModel m;
        m.directory = "/data/sweep";
        m.image_template = "insulin_1_###.img";
        m.frame_range.assign({1, 10}, MetadataSource::CommandLine);
        m.wavelength_A.assign(0.979, MetadataSource::ImageHeader);
        m.distance_mm.assign(160.0, MetadataSource::ImageHeader);
        m.pixel_mm.assign({0.0816, 0.0816}, MetadataSource::ImageHeader);
        m.beam_mm.assign({100.0, 50.0}, MetadataSource::XdsInp);
        m.oscillation_deg.assign(1.0, MetadataSource::ImageHeader);
        m.cell.assign({78.19, 78.19, 78.19, 90.0, 90.0, 90.0}, MetadataSource::XdsInp);
        m.symmetry.assign(197, MetadataSource::XdsInp);
        const auto path = mosflm_test::output_dir("xds_inp_write") / "XDS.INP";
        write_xds_inp(path, m);
        const auto g = XdsInpImporter::parse_file(path);
        CHECK(g.has_template);
        CHECK(g.image_template == "insulin_1_###.img");
        CHECK(g.directory.string() == "/data/sweep");
        CHECK(g.has_process_range);
        CHECK_EQ(g.process_first, 1);
        CHECK_EQ(g.process_last, 10);
        CHECK(g.has_wavelength);
        CHECK(std::fabs(g.wavelength_A - 0.979) < 1e-6);
        CHECK(g.has_distance);
        CHECK(std::fabs(g.distance_mm - 160.0) < 1e-6);
        CHECK(g.has_pixel);
        CHECK(g.has_beam);
        CHECK(std::fabs(g.beam_x_mm - 100.0) < 1e-2);
        CHECK(std::fabs(g.beam_y_mm - 50.0) < 1e-2);
        CHECK(g.has_oscillation);
        CHECK(std::fabs(g.oscillation_deg - 1.0) < 1e-6);
        CHECK(g.has_cell);
        CHECK(std::fabs(g.cell_a - 78.19) < 1e-6);
        CHECK(g.has_symmetry);
        CHECK_EQ(g.symmetry_number, 197);
    }

    return mosflm_test::finish("xinfo_xds_import");
}
