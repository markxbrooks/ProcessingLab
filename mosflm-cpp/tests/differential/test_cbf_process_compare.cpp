/**
 * Differential: CBF PROCESS — F77 ``ipmosflm`` vs ProcessingLab
 * ``ProcessPipeline``.
 *
 * When a Pilatus/miniCBF frame and a Fortran oracle binary are available,
 * runs spot-finding + DPS AUTOINDEX + PROCESS under ipmosflm, then integrates
 * the same cell / orientation matrix in C++ and compares unique-HKL sets
 * (ASU-aware recall gate).
 *
 * Soft-skip (exit 0) if either prerequisite is missing so CI without CCP4 /
 * oracle still passes. Override paths with:
 *   MOSFLM_CBF_IMAGE   — path to ``*.cbf``
 *   MOSFLM_IPMOSFLM    — path to F77 ``ipmosflm`` (e.g. ``bin/ipmosflm.oracle``)
 */

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/crystal/MtzAdapter.hpp"
#include "mosflm/diffraction/ProcessPipeline.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/image/DiffractionImage.hpp"
#include "test_util.hpp"

#include <gemmi/symmetry.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace fs = std::filesystem;

using mosflm::OrientationMatrix;
using mosflm::crystal::CrystalMetadata;
using mosflm::crystal::MtzAdapter;
using mosflm::diffraction::ProcessPipeline;
using mosflm::diffraction::ProcessRequest;
using mosflm::image::CbfImageReader;

namespace {

using Hkl = std::tuple<int, int, int>;

struct Newmat {
    Eigen::Matrix3d amat{Eigen::Matrix3d::Identity()};
    double a{}, b{}, c{}, alpha{}, beta{}, gamma{};
};

void skip_pass(const std::string& why) {
    std::cout << "cbf_process_compare: SKIP (" << why << ")\n";
    std::exit(0);
}

int run_cmd(const std::string& cmd) {
    const int rc = std::system(cmd.c_str());
    if (rc == -1) {
        return 127;
    }
    if ((rc & 0xff) == 0) {
        return (rc >> 8) & 0xff;
    }
    return rc;
}

std::string slurp(const fs::path& path) {
    std::ifstream in(path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

fs::path find_cbf() {
    if (const char* env = std::getenv("MOSFLM_CBF_IMAGE")) {
        return env;
    }
    const char* candidates[] = {
        "/Users/brooks/projects/mxflask/mxpipe/test_data/"
        "KIT-02AT05_w1_1_0001.cbf",
        "../mxflask/mxpipe/test_data/KIT-02AT05_w1_1_0001.cbf",
        "../../mxflask/mxpipe/test_data/KIT-02AT05_w1_1_0001.cbf",
    };
    for (const char* c : candidates) {
        if (fs::exists(c)) {
            return c;
        }
    }
    return {};
}

fs::path find_ipmosflm() {
    if (const char* env = std::getenv("MOSFLM_IPMOSFLM")) {
        return env;
    }
    const char* candidates[] = {
        "../bin/ipmosflm.oracle",
        "../bin/ipmosflm",
        "bin/ipmosflm.oracle",
        "bin/ipmosflm",
        "/Applications/ccp4-9/bin/ipmosflm",
    };
    for (const char* c : candidates) {
        if (fs::exists(c)) {
            return c;
        }
    }
    if (const char* ccp4 = std::getenv("CCP4")) {
        const auto p = fs::path(ccp4) / "bin/ipmosflm";
        if (fs::exists(p)) {
            return p;
        }
    }
    return {};
}

std::string dyld_prefix() {
    // Local oracle links CCP4 + Homebrew HDF5; keep existing DYLD if set.
    const char* existing = std::getenv("DYLD_LIBRARY_PATH");
    std::string path =
        "/Applications/ccp4-9/lib:/opt/homebrew/opt/hdf5/lib";
    if (existing != nullptr && existing[0] != '\0') {
        path = std::string(existing) + ":" + path;
    }
    return "DYLD_LIBRARY_PATH=" + path + " ";
}

Newmat read_newmat(const fs::path& path) {
    std::ifstream in(path);
    CHECK_MSG(static_cast<bool>(in), "cannot open NEWMAT");
    Newmat out;
    double v[9];
    for (double& x : v) {
        CHECK_MSG(static_cast<bool>(in >> x), "NEWMAT matrix truncated");
    }
    out.amat << v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8];
    // Skip misset line + U matrix (3+3 lines of 3 floats), then cell.
    double discard = 0.0;
    for (int i = 0; i < 3 + 9; ++i) {
        CHECK_MSG(static_cast<bool>(in >> discard), "NEWMAT mid section");
    }
    CHECK_MSG(
        static_cast<bool>(
            in >> out.a >> out.b >> out.c >> out.alpha >> out.beta >> out.gamma
        ),
        "NEWMAT cell truncated"
    );
    return out;
}

std::set<Hkl> unique_hkl_mtz(const fs::path& mtz) {
    std::set<Hkl> out;
    for (const auto& r : MtzAdapter::read_simple(mtz)) {
        out.emplace(r.h, r.k, r.l);
    }
    return out;
}

std::set<Hkl> unique_hkl_simple(
    const std::vector<mosflm::crystal::SimpleReflection>& refs
) {
    std::set<Hkl> out;
    for (const auto& r : refs) {
        out.emplace(r.h, r.k, r.l);
    }
    return out;
}

std::set<Hkl> expand_symmetry(
    const std::set<Hkl>& in,
    const gemmi::SpaceGroup* sg
) {
    std::set<Hkl> out;
    const auto ops = sg->operations();
    for (const auto& hkl : in) {
        const gemmi::Op::Miller m{
            std::get<0>(hkl), std::get<1>(hkl), std::get<2>(hkl)
        };
        for (const auto& op : ops) {
            const gemmi::Op::Miller m2 = op.apply_to_hkl(m);
            out.emplace(m2[0], m2[1], m2[2]);
            out.emplace(-m2[0], -m2[1], -m2[2]);
        }
    }
    return out;
}

std::size_t covered_by(
    const std::set<Hkl>& query,
    const std::set<Hkl>& expanded
) {
    std::size_t n = 0;
    for (const auto& h : query) {
        if (expanded.count(h) != 0U) {
            ++n;
        }
    }
    return n;
}

}  // namespace

int main() {
    const fs::path cbf_raw = find_cbf();
    if (cbf_raw.empty() || !fs::exists(cbf_raw)) {
        skip_pass("no CBF image (set MOSFLM_CBF_IMAGE)");
    }
    const fs::path oracle_raw = find_ipmosflm();
    if (oracle_raw.empty() || !fs::exists(oracle_raw)) {
        skip_pass("no F77 ipmosflm (set MOSFLM_IPMOSFLM or build fortran-oracle)");
    }
    // Absolute paths: subprocesses ``cd`` into a temp work directory.
    const fs::path cbf = fs::weakly_canonical(cbf_raw);
    const fs::path oracle = fs::weakly_canonical(oracle_raw);

    const auto image = CbfImageReader{}.read(cbf);
    CHECK_MSG(image.metadata.width > 100, "CBF width too small");
    CHECK_MSG(image.pixels.size() > 0, "CBF has no pixels");

    // MOSFLM Pilatus header beam is (slow, fast) vs Beam_xy (fast, slow).
    const double beam_x_mm = image.metadata.beam_y_mm > 0.0
                                 ? image.metadata.beam_y_mm
                                 : image.metadata.beam_x_mm;
    const double beam_y_mm = image.metadata.beam_x_mm > 0.0
                                 ? image.metadata.beam_x_mm
                                 : image.metadata.beam_y_mm;
    const double dist = image.metadata.detector_distance_mm > 0.0
                            ? image.metadata.detector_distance_mm
                            : 260.0;
    const double wave = image.metadata.wavelength_angstrom > 0.0
                            ? image.metadata.wavelength_angstrom
                            : 0.98;
    const double px = image.metadata.pixel_size_x_mm > 0.0
                          ? image.metadata.pixel_size_x_mm
                          : 0.172;
    const double py = image.metadata.pixel_size_y_mm > 0.0
                          ? image.metadata.pixel_size_y_mm
                          : 0.172;
    const double phi0 = image.metadata.oscillation_start_deg;
    const double dphi = image.metadata.oscillation_range_deg > 0.0
                            ? image.metadata.oscillation_range_deg
                            : 0.1;

    const fs::path work = fs::temp_directory_path() / "mosflm_cbf_process_cmp";
    fs::remove_all(work);
    fs::create_directories(work);
    const fs::path img_link = work / "frame_001.cbf";
    fs::create_symlink(fs::absolute(cbf), img_link);

    // Phase 1: spot find → frame.spt
    {
        std::ofstream in(work / "phase1.in");
        in << "DETECTOR PILATUS\n"
           << "DIRECTORY " << work.string() << "\n"
           << "TEMPLATE frame_###.cbf\n"
           << "BEAM " << beam_x_mm << " " << beam_y_mm << "\n"
           << "DISTANCE " << dist << "\n"
           << "WAVELENGTH " << wave << "\n"
           << "PIXEL " << px << " " << py << "\n"
           << "IMAGE 1\n"
           << "AUTOINDEX\n"
           << "GO\n"
           << "EXIT\n";
    }
    {
        const std::string cmd = "cd " + work.string() + " && " + dyld_prefix() +
                                oracle.string() + " < phase1.in > phase1.log 2>&1";
        CHECK_EQ(run_cmd(cmd), 0);
    }

    // Spot list name follows TEMPLATE stem ("frame").
    const fs::path spt = work / "frame.spt";
    if (!fs::exists(spt)) {
        const std::string log1 = slurp(work / "phase1.log");
        CHECK_MSG(false, "missing frame.spt; phase1.log tail:\n" +
                             log1.substr(log1.size() > 800 ? log1.size() - 800
                                                           : 0));
        return mosflm_test::finish("cbf_process_compare");
    }
    CHECK_MSG(fs::file_size(spt) > 100, "frame.spt too small");

    // Phase 2: DPS index → NEWMAT
    {
        std::ofstream in(work / "phase2.in");
        in << "DETECTOR PILATUS\n"
           << "DIRECTORY " << work.string() << "\n"
           << "TEMPLATE frame_###.cbf\n"
           << "BEAM " << beam_x_mm << " " << beam_y_mm << "\n"
           << "DISTANCE " << dist << "\n"
           << "WAVELENGTH " << wave << "\n"
           << "PIXEL " << px << " " << py << "\n"
           << "AUTOINDEX DPS THRESHOLD 10 FILE " << spt.filename().string()
           << " IMAGE 1\n"
           << "GO\n"
           << "EXIT\n";
    }
    {
        const std::string cmd = "cd " + work.string() + " && " + dyld_prefix() +
                                oracle.string() + " < phase2.in > phase2.log 2>&1";
        CHECK_EQ(run_cmd(cmd), 0);
    }
    const fs::path newmat = work / "NEWMAT";
    CHECK_MSG(fs::exists(newmat), "AUTOINDEX did not write NEWMAT");
    const Newmat nm = read_newmat(newmat);

    // Phase 3: F77 PROCESS → oracle MTZ
    const fs::path oracle_mtz = work / "oracle_cbf.mtz";
    {
        std::ostringstream mat;
        mat << nm.amat(0, 0) << " " << nm.amat(0, 1) << " " << nm.amat(0, 2)
            << " " << nm.amat(1, 0) << " " << nm.amat(1, 1) << " "
            << nm.amat(1, 2) << " " << nm.amat(2, 0) << " " << nm.amat(2, 1)
            << " " << nm.amat(2, 2);
        std::ofstream in(work / "phase3.in");
        in << "DETECTOR PILATUS\n"
           << "DIRECTORY " << work.string() << "\n"
           << "TEMPLATE frame_###.cbf\n"
           << "BEAM " << beam_x_mm << " " << beam_y_mm << "\n"
           << "DISTANCE " << dist << "\n"
           << "WAVELENGTH " << wave << "\n"
           << "PIXEL " << px << " " << py << "\n"
           << "CELL " << nm.a << " " << nm.b << " " << nm.c << " " << nm.alpha
           << " " << nm.beta << " " << nm.gamma << "\n"
           << "SYMMETRY 1\n"
           << "MOSAIC 0.50\n"
           << "MATRIX " << mat.str() << "\n"
           << "PROCESS 1 TO 1 START " << phi0 << " ANGLE " << dphi << "\n"
           << "GO\n"
           << "EXIT\n";
    }
    {
        const std::string cmd = "cd " + work.string() + " && " + dyld_prefix() +
                                oracle.string() + " HKLOUT " +
                                oracle_mtz.string() +
                                " < phase3.in > phase3.log 2>&1";
        const int rc = run_cmd(cmd);
        CHECK_MSG(rc == 0 && fs::exists(oracle_mtz) &&
                      fs::file_size(oracle_mtz) > 0,
                  "F77 PROCESS failed (see phase3.log in work dir)");
    }

    // C++ PROCESS with the same cell / AMAT / beam.
    mosflm::DetectorGeometry det = image.metadata.to_detector_geometry();
    det.beam_x_mm = 0.0;
    det.beam_y_mm = 0.0;
    if (det.distance_mm <= 0.0) {
        det.distance_mm = dist;
    }
    if (det.wavelength_A <= 0.0) {
        det.wavelength_A = wave;
    }

    mosflm::diffraction::ScanRange scan{phi0, phi0 + dphi};
    mosflm::diffraction::PredictionSettings prediction;
    prediction.dmin_A = 1.5;
    prediction.ewald_tolerance_Ainv = 0.0025;
    prediction.mosaic_spread_deg = 0.50;
    prediction.phi_step_deg = 0.025;

    mosflm::diffraction::IntegrationSettings integration;
    integration.box_half_px = 3;
    integration.bg_frame_px = 2;
    integration.apply_lp = true;
    integration.mosflm_detector_axes = true;
    integration.image_beam_x_mm = beam_x_mm;
    integration.image_beam_y_mm = beam_y_mm;

    mosflm::diffraction::MeasurementBoxSettings measurement;
    measurement.separation_x_mm = 0.7;
    measurement.separation_y_mm = 1.0;
    measurement.bg_rim_px = 2;
    measurement.mark_overlaps = true;
    measurement.integrate_overlaps = true;

    mosflm::diffraction::StandardProfileSettings profile;
    profile.box = measurement;
    profile.learn_min_isigi = 5.0;
    profile.min_learners = 8;
    profile.prefer_profile = true;
    profile.n_bins_x = 5;
    profile.n_bins_y = 5;

    const ProcessRequest req{
        CrystalMetadata::from_parameters(
            nm.a, nm.b, nm.c, nm.alpha, nm.beta, nm.gamma, "P 1"
        ),
        OrientationMatrix::from_mosflm_orientation_matrix(nm.amat),
        det,
        scan,
        prediction,
        integration,
        measurement,
        profile,
        true,
        true,
        false,
        {},
        0.0,
        image
    };
    const auto cpp = ProcessPipeline{}.run(req);
    CHECK_MSG(cpp.predictions.size() > 10, "C++ predicted too few spots");
    CHECK_MSG(
        cpp.mtz_reflections.size() > 5,
        "C++ integrated too few reflections"
    );

    const auto oracle_hkl = unique_hkl_mtz(oracle_mtz);
    const auto cpp_hkl = unique_hkl_simple(cpp.mtz_reflections);
    CHECK_MSG(oracle_hkl.size() > 10, "oracle HKL set too small");

    const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_name("P 1");
    CHECK_MSG(sg != nullptr, "missing P1");
    const auto cpp_expanded = expand_symmetry(cpp_hkl, sg);
    const std::size_t asu_hits = covered_by(oracle_hkl, cpp_expanded);
    const double asu_recall =
        static_cast<double>(asu_hits) /
        static_cast<double>(oracle_hkl.size());

    std::size_t exact = 0;
    for (const auto& h : oracle_hkl) {
        if (cpp_hkl.count(h) != 0U) {
            ++exact;
        }
    }
    const double exact_recall =
        static_cast<double>(exact) / static_cast<double>(oracle_hkl.size());

    std::cout << "cbf_process_compare: cbf=" << cbf
              << " oracle_bin=" << oracle << "\n"
              << "  cell=" << nm.a << " " << nm.b << " " << nm.c << " "
              << nm.alpha << " " << nm.beta << " " << nm.gamma << "\n"
              << "  f77_unique=" << oracle_hkl.size()
              << " cpp_unique=" << cpp_hkl.size()
              << " cpp_preds=" << cpp.predictions.size()
              << " exact_recall=" << exact_recall
              << " asu_recall=" << asu_recall << '\n';

    // Soft gate: CBF path is newer than insulin SMV; require meaningful overlap.
    CHECK_MSG(asu_recall >= 0.40, "ASU-aware HKL recall below 0.40 vs F77 CBF");

    fs::remove_all(work);
    return mosflm_test::finish("cbf_process_compare");
}
