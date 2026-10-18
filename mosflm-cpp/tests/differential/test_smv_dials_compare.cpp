/**
 * Differential: insulin SMV PROCESS — ProcessingLab (C++) vs DIALS.
 *
 * Runs ``cpp_process.in`` through ProcessingLab, then a DIALS pipeline
 * (import → find_spots → index with known I23/cell → integrate → export MTZ)
 * on the same CCP4 insulin SMV frame.
 *
 * DIALS re-indexes, so HKL labels need not match MOSFLM; the primary science
 * gate is detector-plane proximity (XDET/YDET). Cell/SG and soft count /
 * I-ratio gates mirror ``smv_process_compare``.
 *
 * Soft-skip (exit 0) if ProcessingLab, DIALS, or the image is missing.
 * Override paths with:
 *   MOSFLM_PROCESSINGLAB — C++ CONTROL binary (ProcessingLab)
 *   MOSFLM_DIALS_BIN       — directory containing ``dials.import`` etc.
 *   MOSFLM_IMAGE           — path to ``insulin_1_001.img``
 *   CCP4                   — used to locate DIALS under ``$CCP4/bin``
 */

#include "mosflm/crystal/MtzAdapter.hpp"
#include "test_util.hpp"

#include <gemmi/mtz.hpp>
#include <gemmi/symmetry.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

using mosflm::crystal::MtzAdapter;

namespace {

using Hkl = std::tuple<int, int, int>;
using Xy = std::pair<double, double>;

void skip_pass(const std::string& why) {
    std::cout << "smv_dials_compare: SKIP (" << why << ")\n";
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

fs::path find_integrationlab() {
    if (const char* env = std::getenv("MOSFLM_PROCESSINGLAB")) {
        return env;
    }
    if (const char* env = std::getenv("MOSFLM_THIN_IPMOSFLM")) {
        return env;
    }
    const char* rels[] = {
        "build/ProcessingLab",
        "build/ProcessingLab.exe",
        "./ProcessingLab",
        "./ProcessingLab.exe",
        "../bin/ProcessingLab",
        "../bin/ProcessingLab.exe",
        "bin/ProcessingLab",
    };
    for (const char* rel : rels) {
        if (fs::exists(rel)) {
            return rel;
        }
    }
    return {};
}

fs::path find_dials_bin_dir() {
    if (const char* env = std::getenv("MOSFLM_DIALS_BIN")) {
        const fs::path p = env;
        if (fs::exists(p / "dials.import")) {
            return p;
        }
    }
    if (const char* ccp4 = std::getenv("CCP4")) {
        const fs::path p = fs::path(ccp4) / "bin";
        if (fs::exists(p / "dials.import")) {
            return p;
        }
    }
    const char* candidates[] = {
        "/Applications/ccp4-9/bin",
        "/apps/ccp4/ccp4-9/bin",
    };
    for (const char* c : candidates) {
        if (fs::exists(fs::path(c) / "dials.import")) {
            return c;
        }
    }
    return {};
}

fs::path find_insulin_img() {
    if (const char* env = std::getenv("MOSFLM_IMAGE")) {
        return env;
    }
    if (const char* ccp4 = std::getenv("CCP4")) {
        const auto p =
            fs::path(ccp4) / "examples/data/insulin_1_001.img";
        if (fs::exists(p)) {
            return p;
        }
    }
    const char* fallbacks[] = {
        "/Applications/ccp4-9/examples/data/insulin_1_001.img",
        "/apps/ccp4/ccp4-9/examples/data/insulin_1_001.img",
    };
    for (const char* c : fallbacks) {
        if (fs::exists(c)) {
            return c;
        }
    }
    return {};
}

/// CCP4 DIALS prefers a clean PYTHONPATH (venv site-packages can break import).
std::string dials_env_prefix(const fs::path& dials_bin) {
    std::string ccp4 = "/Applications/ccp4-9";
    if (const char* env = std::getenv("CCP4")) {
        if (env[0] != '\0') {
            ccp4 = env;
        }
    }
    return "env -u PYTHONPATH -u PYTHONHOME CCP4=" + ccp4 + " PATH=" +
           dials_bin.string() + ":\"$PATH\" ";
}

std::set<Hkl> unique_hkl_mtz(const fs::path& mtz) {
    std::set<Hkl> out;
    for (const auto& r : MtzAdapter::read_simple(mtz)) {
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

struct Spot {
    double xdet{};
    double ydet{};
    double intensity{};
    double sigma{};
};

std::vector<Spot> read_spots(const fs::path& path) {
    gemmi::Mtz mtz;
    mtz.read_file(path.string());
    const auto* X = mtz.column_with_label("XDET");
    const auto* Y = mtz.column_with_label("YDET");
    const auto* I = mtz.column_with_label("I");
    const auto* S = mtz.column_with_label("SIGI");
    CHECK_MSG(X && Y && I && S, "MTZ missing XDET/YDET/I/SIGI");
    std::vector<Spot> out;
    out.reserve(static_cast<std::size_t>(mtz.nreflections));
    for (int i = 0; i < mtz.nreflections; ++i) {
        Spot s;
        s.xdet = (*X)[static_cast<std::size_t>(i)];
        s.ydet = (*Y)[static_cast<std::size_t>(i)];
        s.intensity = (*I)[static_cast<std::size_t>(i)];
        s.sigma = (*S)[static_cast<std::size_t>(i)];
        out.push_back(s);
    }
    return out;
}

double nearest_dist_px(const Spot& q, const std::vector<Spot>& cloud) {
    double best = std::numeric_limits<double>::infinity();
    for (const auto& s : cloud) {
        const double dx = q.xdet - s.xdet;
        const double dy = q.ydet - s.ydet;
        const double d2 = dx * dx + dy * dy;
        if (d2 < best) {
            best = d2;
        }
    }
    return std::sqrt(best);
}

const Spot* nearest_spot(const Spot& q, const std::vector<Spot>& cloud) {
    const Spot* best = nullptr;
    double best_d2 = std::numeric_limits<double>::infinity();
    for (const auto& s : cloud) {
        const double dx = q.xdet - s.xdet;
        const double dy = q.ydet - s.ydet;
        const double d2 = dx * dx + dy * dy;
        if (d2 < best_d2) {
            best_d2 = d2;
            best = &s;
        }
    }
    return best;
}

}  // namespace

int main() {
    const fs::path cpp_bin_raw = find_integrationlab();
    const fs::path dials_bin_raw = find_dials_bin_dir();
    const fs::path img_raw = find_insulin_img();
    const fs::path cpp_script =
        "testdata/images/known_indexed/cpp_process.in";

    if (cpp_bin_raw.empty() || !fs::exists(cpp_bin_raw)) {
        skip_pass("no ProcessingLab (set MOSFLM_PROCESSINGLAB)");
    }
    if (dials_bin_raw.empty()) {
        skip_pass("no DIALS (set MOSFLM_DIALS_BIN or CCP4)");
    }
    if (img_raw.empty() || !fs::exists(img_raw)) {
        skip_pass("no insulin SMV (set MOSFLM_IMAGE or CCP4)");
    }
    CHECK_MSG(fs::exists(cpp_script), "missing cpp_process.in");

    const fs::path cpp_bin = fs::weakly_canonical(cpp_bin_raw);
    const fs::path dials_bin = fs::weakly_canonical(dials_bin_raw);
    const fs::path img = fs::weakly_canonical(img_raw);

    const fs::path work =
        fs::temp_directory_path() / "mosflm_smv_dials_cmp";
    fs::remove_all(work);
    fs::create_directories(work);

    // --- ProcessingLab PROCESS ------------------------------------------------
    const fs::path cpp_mtz = work / "cpp.mtz";
    {
        const std::string cmd =
            cpp_bin.string() + " HKLOUT " + cpp_mtz.string() + " --image " +
            img.string() + " " + cpp_script.string() + " >" +
            (work / "cpp.log").string() + " 2>&1";
        CHECK_EQ(run_cmd(cmd), 0);
        CHECK_MSG(
            fs::exists(cpp_mtz) && fs::file_size(cpp_mtz) > 0,
            "C++ PROCESS produced no MTZ"
        );
    }

    // --- DIALS: import → find_spots → index → integrate → export -------------
    const fs::path dials_work = work / "dials";
    fs::create_directories(dials_work);
    const fs::path dials_mtz = dials_work / "dials.mtz";
    const std::string env = dials_env_prefix(dials_bin);
    {
        const std::string cd = "cd " + dials_work.string() + " && " + env;
        CHECK_EQ(
            run_cmd(
                cd + "dials.import " + img.string() + " > import.log 2>&1"
            ),
            0
        );
        CHECK_EQ(
            run_cmd(cd + "dials.find_spots imported.expt nproc=2 > find.log 2>&1"
            ),
            0
        );
        CHECK_EQ(
            run_cmd(
                cd +
                "dials.index imported.expt strong.refl "
                "unit_cell=78.19,78.19,78.19,90,90,90 space_group=I23 "
                "nproc=2 > index.log 2>&1"
            ),
            0
        );
        CHECK_EQ(
            run_cmd(
                cd +
                "dials.integrate indexed.expt indexed.refl nproc=2 "
                "> integrate.log 2>&1"
            ),
            0
        );
        CHECK_EQ(
            run_cmd(
                cd +
                "dials.export integrated.expt integrated.refl format=mtz "
                "mtz.hklout=dials.mtz > export.log 2>&1"
            ),
            0
        );
        CHECK_MSG(
            fs::exists(dials_mtz) && fs::file_size(dials_mtz) > 0,
            "DIALS export produced no dials.mtz (see dials/*.log)"
        );
    }

    gemmi::Mtz cpp_hdr;
    gemmi::Mtz dials_hdr;
    cpp_hdr.read_file(cpp_mtz.string());
    dials_hdr.read_file(dials_mtz.string());

    CHECK_MSG(cpp_hdr.nreflections > 50, "C++ MTZ too small");
    CHECK_MSG(dials_hdr.nreflections > 50, "DIALS MTZ too small");
    // DIALS refines cell slightly; allow ~1 Å / 1°.
    CHECK_CLOSE(cpp_hdr.cell.a, dials_hdr.cell.a, 1.0);
    CHECK_CLOSE(cpp_hdr.cell.b, dials_hdr.cell.b, 1.0);
    CHECK_CLOSE(cpp_hdr.cell.c, dials_hdr.cell.c, 1.0);
    CHECK_CLOSE(cpp_hdr.cell.alpha, dials_hdr.cell.alpha, 1.0);
    CHECK_CLOSE(cpp_hdr.cell.beta, dials_hdr.cell.beta, 1.0);
    CHECK_CLOSE(cpp_hdr.cell.gamma, dials_hdr.cell.gamma, 1.0);
    CHECK_MSG(cpp_hdr.spacegroup != nullptr, "C++ MTZ missing space group");
    CHECK_MSG(dials_hdr.spacegroup != nullptr, "DIALS MTZ missing space group");
    CHECK_EQ(cpp_hdr.spacegroup->number, dials_hdr.spacegroup->number);

    CHECK_MSG(
        cpp_hdr.column_with_label("XDET") != nullptr &&
            cpp_hdr.column_with_label("YDET") != nullptr,
        "C++ MTZ missing XDET/YDET"
    );
    CHECK_MSG(
        dials_hdr.column_with_label("XDET") != nullptr &&
            dials_hdr.column_with_label("YDET") != nullptr,
        "DIALS MTZ missing XDET/YDET"
    );

    const auto cpp_hkl = unique_hkl_mtz(cpp_mtz);
    const auto dials_hkl = unique_hkl_mtz(dials_mtz);
    const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_number(197);
    CHECK_MSG(sg != nullptr, "missing I23");
    const auto cpp_expanded = expand_symmetry(cpp_hkl, sg);
    const std::size_t asu_hits = covered_by(dials_hkl, cpp_expanded);
    const double asu_recall =
        dials_hkl.empty()
            ? 0.0
            : static_cast<double>(asu_hits) /
                  static_cast<double>(dials_hkl.size());

    const double count_ratio =
        static_cast<double>(cpp_hdr.nreflections) /
        static_cast<double>(dials_hdr.nreflections);

    // Detector-plane proximity (robust when DIALS UB ≠ MOSFLM MATRIX).
    const auto cpp_spots = read_spots(cpp_mtz);
    const auto dials_spots = read_spots(dials_mtz);
    std::vector<double> nearest;
    nearest.reserve(dials_spots.size());
    std::size_t within_30 = 0;
    std::vector<double> i_ratios;
    for (const auto& d : dials_spots) {
        const double dist = nearest_dist_px(d, cpp_spots);
        nearest.push_back(dist);
        if (dist <= 30.0) {
            ++within_30;
        }
        if (d.intensity > 0.0 && d.sigma > 0.0) {
            const Spot* c = nearest_spot(d, cpp_spots);
            if (c != nullptr && dist <= 30.0 && c->intensity > 0.0 &&
                c->sigma > 0.0) {
                i_ratios.push_back(c->intensity / d.intensity);
            }
        }
    }
    std::sort(nearest.begin(), nearest.end());
    const double med_xy =
        nearest.empty() ? 0.0 : nearest[nearest.size() / 2];
    const double frac_xy30 =
        dials_spots.empty()
            ? 0.0
            : static_cast<double>(within_30) /
                  static_cast<double>(dials_spots.size());
    double med_i_ratio = 0.0;
    if (!i_ratios.empty()) {
        std::sort(i_ratios.begin(), i_ratios.end());
        med_i_ratio = i_ratios[i_ratios.size() / 2];
    }

    std::cout << "smv_dials_compare: C++=" << cpp_bin << "\n"
              << "smv_dials_compare: DIALS_BIN=" << dials_bin << "\n"
              << "smv_dials_compare: image=" << img << "\n"
              << "  cell_cpp=" << cpp_hdr.cell.a << " cell_dials="
              << dials_hdr.cell.a << " sg=" << cpp_hdr.spacegroup->number
              << "\n"
              << "  dials_nref=" << dials_hdr.nreflections
              << " dials_unique=" << dials_hkl.size()
              << " cpp_nref=" << cpp_hdr.nreflections
              << " cpp_unique=" << cpp_hkl.size()
              << " count_ratio=" << count_ratio
              << " asu_recall=" << asu_recall
              << " frac_xy<=30px=" << frac_xy30
              << " med_xy_px=" << med_xy
              << " matched_I=" << i_ratios.size()
              << " med_I_ratio=" << med_i_ratio << '\n';

    CHECK_MSG(
        count_ratio >= 0.5 && count_ratio <= 5.0,
        "C++/DIALS row-count ratio outside soft band [0.5, 5]"
    );
    // Indexing differs → HKL ASU recall is informational / soft floor only.
    CHECK_MSG(
        asu_recall >= 0.05,
        "ASU HKL recall of DIALS indices below soft floor 0.05"
    );
    CHECK_MSG(
        frac_xy30 >= 0.40,
        "fewer than 40% of DIALS spots within 30 px of a C++ XDET/YDET"
    );
    CHECK_MSG(med_xy < 40.0, "median nearest XDET distance above 40 px");
    if (i_ratios.size() >= 50) {
        CHECK_MSG(
            med_i_ratio > 0.05 && med_i_ratio < 50.0,
            "median XDET-matched I ratio out of soft band"
        );
    }

    fs::remove_all(work);
    return mosflm_test::finish("smv_dials_compare");
}
