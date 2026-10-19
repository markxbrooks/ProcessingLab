/**
 * Differential: insulin SMV PROCESS — ProcessingLab (C++) vs F77 ``ipmosflm``.
 *
 * Runs the known-indexed CONTROL scripts (``cpp_process.in`` /
 * ``ipmosflm_process.in``) against the CCP4 insulin frame, then compares the
 * resulting MTZs the same way as a manual mtzdmp/gemmi check:
 *   - unit cell / space group agreement
 *   - unique-HKL exact overlap + ASU-aware recall
 *   - soft row-count ratio and median I-ratio on ASU-matched pairs
 *
 * Soft-skip (exit 0) if either binary or the image is missing so CI without
 * CCP4 / oracle still passes. Override paths with:
 *   MOSFLM_PROCESSINGLAB — C++ CONTROL binary (ProcessingLab)
 *   MOSFLM_IPMOSFLM       — F77 ``ipmosflm`` / ``ipmosflm.oracle``
 *   MOSFLM_IMAGE          — path to ``insulin_1_001.img``
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
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace fs = std::filesystem;

using mosflm::crystal::MtzAdapter;

namespace {

using Hkl = std::tuple<int, int, int>;

void skip_pass(const std::string& why) {
    std::cout << "smv_process_compare: SKIP (" << why << ")\n";
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
        "../bin/IntegrationLab",
        "bin/IntegrationLab",
    };
    for (const char* rel : rels) {
        if (fs::exists(rel)) {
            return rel;
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
        "/apps/ccp4/ccp4-9/bin/ipmosflm",
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

std::string dyld_prefix() {
#if defined(__APPLE__)
    const char* existing = std::getenv("DYLD_LIBRARY_PATH");
    std::string path =
        "/Applications/ccp4-9/lib:/opt/homebrew/opt/hdf5/lib";
    if (const char* ccp4 = std::getenv("CCP4")) {
        path = (fs::path(ccp4) / "lib").string() + ":" + path;
    }
    if (existing != nullptr && existing[0] != '\0') {
        path = std::string(existing) + ":" + path;
    }
    return "DYLD_LIBRARY_PATH=" + path + " ";
#else
    return {};
#endif
}

std::set<Hkl> unique_hkl_mtz(const fs::path& mtz) {
    std::set<Hkl> out;
    for (const auto& r : MtzAdapter::read_simple(mtz)) {
        out.emplace(r.h, r.k, r.l);
    }
    return out;
}

std::map<Hkl, std::pair<double, double>> intensity_map(const fs::path& mtz) {
    std::map<Hkl, std::pair<double, double>> out;
    for (const auto& r : MtzAdapter::read_simple(mtz)) {
        const Hkl key{r.h, r.k, r.l};
        const auto it = out.find(key);
        if (it == out.end() || r.intensity > it->second.first) {
            out[key] = {r.intensity, r.sigma};
        }
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

std::map<Hkl, std::pair<double, double>> expand_intensity(
    const std::map<Hkl, std::pair<double, double>>& in,
    const gemmi::SpaceGroup* sg
) {
    std::map<Hkl, std::pair<double, double>> out;
    const auto ops = sg->operations();
    for (const auto& [hkl, spot] : in) {
        const gemmi::Op::Miller m{
            std::get<0>(hkl), std::get<1>(hkl), std::get<2>(hkl)
        };
        for (const auto& op : ops) {
            const gemmi::Op::Miller m2 = op.apply_to_hkl(m);
            for (const Hkl cand :
                 {Hkl{m2[0], m2[1], m2[2]}, Hkl{-m2[0], -m2[1], -m2[2]}}) {
                const auto it = out.find(cand);
                if (it == out.end() || spot.first > it->second.first) {
                    out[cand] = spot;
                }
            }
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

bool has_column(const gemmi::Mtz& mtz, const char* label) {
    return mtz.column_with_label(label) != nullptr;
}

}  // namespace

int main() {
    const fs::path cpp_bin_raw = find_integrationlab();
    const fs::path f77_bin_raw = find_ipmosflm();
    const fs::path img_raw = find_insulin_img();
    const fs::path cpp_script =
        "testdata/images/known_indexed/cpp_process.in";
    const fs::path f77_script =
        "testdata/images/known_indexed/ipmosflm_process.in";
    const fs::path mat =
        "testdata/images/known_indexed/insulin_1.mat";

    if (cpp_bin_raw.empty() || !fs::exists(cpp_bin_raw)) {
        skip_pass("no ProcessingLab (set MOSFLM_PROCESSINGLAB)");
    }
    if (f77_bin_raw.empty() || !fs::exists(f77_bin_raw)) {
        skip_pass("no F77 ipmosflm (set MOSFLM_IPMOSFLM or build fortran-oracle)");
    }
    if (img_raw.empty() || !fs::exists(img_raw)) {
        skip_pass("no insulin SMV (set MOSFLM_IMAGE or CCP4)");
    }
    CHECK_MSG(fs::exists(cpp_script), "missing cpp_process.in");
    CHECK_MSG(fs::exists(f77_script), "missing ipmosflm_process.in");
    CHECK_MSG(fs::exists(mat), "missing insulin_1.mat");

    const fs::path cpp_bin = fs::weakly_canonical(cpp_bin_raw);
    const fs::path f77_bin = fs::weakly_canonical(f77_bin_raw);
    const fs::path img = fs::weakly_canonical(img_raw);

    const fs::path work =
        fs::temp_directory_path() / "mosflm_smv_process_cmp";
    fs::remove_all(work);
    fs::create_directories(work);

    const fs::path cpp_mtz = work / "cpp.mtz";
    std::size_t cpp_preds_logged = 0;
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
        // Parse "preds=N" from ProcessingLab stdout if present.
        std::ifstream log(work / "cpp.log");
        std::string line;
        while (std::getline(log, line)) {
            const auto pos = line.find("preds=");
            if (pos == std::string::npos) {
                continue;
            }
            cpp_preds_logged = static_cast<std::size_t>(
                std::strtoul(line.c_str() + pos + 6, nullptr, 10)
            );
        }
    }

    const fs::path f77_work = work / "f77";
    fs::create_directories(f77_work);
    fs::copy_file(img, f77_work / "insulin_1_001.img");
    fs::copy_file(mat, f77_work / "insulin_1.mat");
    fs::copy_file(f77_script, f77_work / "ipmosflm_process.in");
    const fs::path f77_mtz = f77_work / "PERF.mtz";
    {
        const std::string cmd =
            "cd " + f77_work.string() + " && " + dyld_prefix() +
            f77_bin.string() + " < ipmosflm_process.in > f77.log 2>&1";
        const int rc = run_cmd(cmd);
        CHECK_MSG(
            rc == 0 && fs::exists(f77_mtz) && fs::file_size(f77_mtz) > 0,
            "F77 PROCESS failed (see f77.log in work dir)"
        );
    }

    gemmi::Mtz cpp_hdr;
    gemmi::Mtz f77_hdr;
    cpp_hdr.read_file(cpp_mtz.string());
    f77_hdr.read_file(f77_mtz.string());

    CHECK_MSG(cpp_hdr.nreflections > 50, "C++ MTZ too small");
    CHECK_MSG(f77_hdr.nreflections > 50, "F77 MTZ too small");
    CHECK_CLOSE(cpp_hdr.cell.a, f77_hdr.cell.a, 1e-2);
    CHECK_CLOSE(cpp_hdr.cell.b, f77_hdr.cell.b, 1e-2);
    CHECK_CLOSE(cpp_hdr.cell.c, f77_hdr.cell.c, 1e-2);
    CHECK_CLOSE(cpp_hdr.cell.alpha, f77_hdr.cell.alpha, 1e-2);
    CHECK_CLOSE(cpp_hdr.cell.beta, f77_hdr.cell.beta, 1e-2);
    CHECK_CLOSE(cpp_hdr.cell.gamma, f77_hdr.cell.gamma, 1e-2);
    CHECK_MSG(cpp_hdr.spacegroup != nullptr, "C++ MTZ missing space group");
    CHECK_MSG(f77_hdr.spacegroup != nullptr, "F77 MTZ missing space group");
    CHECK_EQ(cpp_hdr.spacegroup->number, f77_hdr.spacegroup->number);

    // Full MOSFLM schema when C++ writer is enabled (Phase 3).
    const bool cpp_full_schema =
        has_column(cpp_hdr, "M/ISYM") && has_column(cpp_hdr, "XDET") &&
        has_column(cpp_hdr, "FLAG") && has_column(cpp_hdr, "FRACTIONCALC");
    if (cpp_full_schema) {
        CHECK_MSG(has_column(cpp_hdr, "BATCH"), "missing BATCH");
        CHECK_MSG(has_column(cpp_hdr, "IPR"), "missing IPR");
        CHECK_MSG(has_column(cpp_hdr, "YDET"), "missing YDET");
        CHECK_MSG(has_column(cpp_hdr, "LP"), "missing LP");
        CHECK_MSG(has_column(cpp_hdr, "MPART"), "missing MPART");
        CHECK_MSG(has_column(cpp_hdr, "BGPKRATIOS"), "missing BGPKRATIOS");
    }

    const auto cpp_hkl = unique_hkl_mtz(cpp_mtz);
    const auto f77_hkl = unique_hkl_mtz(f77_mtz);
    CHECK_MSG(cpp_hkl.size() > 50, "C++ unique HKL set too small");
    CHECK_MSG(f77_hkl.size() > 50, "F77 unique HKL set too small");

    const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_number(197);
    CHECK_MSG(sg != nullptr, "missing I23 space group");
    const auto cpp_expanded = expand_symmetry(cpp_hkl, sg);
    const std::size_t asu_hits = covered_by(f77_hkl, cpp_expanded);
    const double asu_recall =
        static_cast<double>(asu_hits) /
        static_cast<double>(f77_hkl.size());

    std::size_t exact = 0;
    for (const auto& h : f77_hkl) {
        if (cpp_hkl.count(h) != 0U) {
            ++exact;
        }
    }
    const double exact_recall =
        static_cast<double>(exact) / static_cast<double>(f77_hkl.size());

    const double count_ratio =
        static_cast<double>(cpp_hdr.nreflections) /
        static_cast<double>(f77_hdr.nreflections);

    // Soft I-ratio on ASU-matched positive pairs.
    const auto f77_I = expand_intensity(intensity_map(f77_mtz), sg);
    const auto cpp_I = intensity_map(cpp_mtz);
    std::vector<double> ratios;
    std::vector<double> abs_z;
    for (const auto& [hkl, spot] : cpp_I) {
        if (spot.second <= 0.0 || spot.first <= 0.0) {
            continue;
        }
        const auto it = f77_I.find(hkl);
        if (it == f77_I.end() || it->second.first <= 0.0 ||
            it->second.second <= 0.0) {
            continue;
        }
        ratios.push_back(spot.first / it->second.first);
        abs_z.push_back(
            std::fabs(spot.first - it->second.first) / it->second.second
        );
    }
    double med_ratio = 0.0;
    double med_abs_z = 0.0;
    if (!ratios.empty()) {
        std::sort(ratios.begin(), ratios.end());
        std::sort(abs_z.begin(), abs_z.end());
        med_ratio = ratios[ratios.size() / 2];
        med_abs_z = abs_z[abs_z.size() / 2];
    }

    std::cout << "smv_process_compare: C++=" << cpp_bin << "\n"
              << "smv_process_compare: F77=" << f77_bin << "\n"
              << "smv_process_compare: image=" << img << "\n"
              << "  cell=" << cpp_hdr.cell.a << " " << cpp_hdr.cell.b << " "
              << cpp_hdr.cell.c << " sg=" << cpp_hdr.spacegroup->number
              << " full_schema=" << (cpp_full_schema ? 1 : 0) << "\n"
              << "  cpp_preds=" << cpp_preds_logged
              << " f77_nref=" << f77_hdr.nreflections
              << " f77_unique=" << f77_hkl.size()
              << " cpp_nref=" << cpp_hdr.nreflections
              << " cpp_unique=" << cpp_hkl.size()
              << " count_ratio=" << count_ratio
              << " exact_overlap=" << exact
              << " exact_recall=" << exact_recall
              << " asu_recall=" << asu_recall
              << " matched_I=" << ratios.size()
              << " med_I_ratio=" << med_ratio
              << " med_|z|=" << med_abs_z << '\n';

    CHECK_MSG(exact > 0, "no exact HKL overlap between C++ and F77 MTZs");
    CHECK_MSG(
        asu_recall >= 0.85,
        "ASU-aware HKL recall of F77 indices below 0.85 vs C++"
    );

    // Soft count gate (target band 0.85–1.15 once REEK funnel matches).
    CHECK_MSG(
        count_ratio >= 0.5 && count_ratio <= 3.0,
        "C++/F77 row-count ratio outside soft band [0.5, 3]"
    );

    if (ratios.size() >= 50) {
        CHECK_MSG(
            med_ratio > 0.3 && med_ratio < 5.0,
            "median ASU-matched I ratio out of band [0.3, 5]"
        );
    }

    // After ASU packing (Phase 4), exact recall should rise.
    if (cpp_full_schema && count_ratio >= 0.85 && count_ratio <= 1.15) {
        CHECK_MSG(
            exact_recall >= 0.80,
            "exact HKL recall below 0.80 after count parity"
        );
        CHECK_MSG(asu_recall >= 0.95, "ASU recall below 0.95 after count parity");
    }

    fs::remove_all(work);
    return mosflm_test::finish("smv_process_compare");
}
