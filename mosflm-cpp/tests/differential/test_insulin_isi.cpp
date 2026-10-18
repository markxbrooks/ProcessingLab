/**
 * Phase D / D2: soft I/SIGI differential — C++ MEAS PROCESS vs Fortran oracle.
 *
 * Matches on space-group / Friedel equivalents. Gates are intentionally soft:
 * XDET residuals and profile fitting (D3) still dominate intensity scale.
 */

#include "mosflm/diffraction/ProcessPipeline.hpp"
#include "mosflm/image/DiffractionImage.hpp"
#include "test_util.hpp"

#include <gemmi/mtz.hpp>
#include <gemmi/symmetry.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <map>
#include <numeric>
#include <tuple>
#include <vector>

using mosflm::diffraction::ProcessPipeline;
using mosflm::image::SmvImageReader;

namespace {

using Hkl = std::tuple<int, int, int>;

std::filesystem::path find_insulin_img() {
    if (const char* env = std::getenv("MOSFLM_IMAGE")) {
        return env;
    }
    if (const char* ccp4 = std::getenv("CCP4")) {
        const auto p =
            std::filesystem::path(ccp4) / "examples/data/insulin_1_001.img";
        if (std::filesystem::exists(p)) {
            return p;
        }
    }
    return "/Applications/ccp4-9/examples/data/insulin_1_001.img";
}

struct OracleSpot {
    double intensity{};
    double sigma{};
};

std::map<Hkl, OracleSpot> read_oracle_intensities(
    const std::filesystem::path& path
) {
    gemmi::Mtz mtz;
    mtz.read_file(path.string());
    const gemmi::Mtz::Column* col_h = mtz.column_with_label("H");
    const gemmi::Mtz::Column* col_k = mtz.column_with_label("K");
    const gemmi::Mtz::Column* col_l = mtz.column_with_label("L");
    const gemmi::Mtz::Column* col_i = mtz.column_with_label("I");
    const gemmi::Mtz::Column* col_s = mtz.column_with_label("SIGI");
    CHECK_MSG(col_h && col_k && col_l && col_i && col_s, "oracle MTZ columns");
    std::map<Hkl, OracleSpot> out;
    for (int i = 0; i < mtz.nreflections; ++i) {
        const Hkl key{
            static_cast<int>(std::lround((*col_h)[static_cast<std::size_t>(i)])),
            static_cast<int>(std::lround((*col_k)[static_cast<std::size_t>(i)])),
            static_cast<int>(std::lround((*col_l)[static_cast<std::size_t>(i)]))
        };
        OracleSpot s{
            (*col_i)[static_cast<std::size_t>(i)],
            (*col_s)[static_cast<std::size_t>(i)]
        };
        const auto it = out.find(key);
        if (it == out.end() || s.intensity > it->second.intensity) {
            out[key] = s;
        }
    }
    return out;
}

std::map<Hkl, OracleSpot> expand_symmetry(
    const std::map<Hkl, OracleSpot>& in,
    const gemmi::SpaceGroup* sg
) {
    std::map<Hkl, OracleSpot> out;
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
                if (it == out.end() || spot.intensity > it->second.intensity) {
                    out[cand] = spot;
                }
            }
        }
    }
    return out;
}

}  // namespace

int main() {
    const std::filesystem::path oracle_path =
        "testdata/oracle/known_indexed/insulin_1_001.mtz";
    CHECK_MSG(std::filesystem::exists(oracle_path), "missing oracle MTZ");
    const auto img_path = find_insulin_img();
    CHECK_MSG(
        std::filesystem::exists(img_path),
        "missing insulin SMV (set MOSFLM_IMAGE or CCP4)"
    );

    const auto image = SmvImageReader{}.read(img_path);
    const auto req = ProcessPipeline::known_indexed_insulin(image);
    CHECK(req.use_meas);
    const auto result = ProcessPipeline{}.run(req);
    CHECK_MSG(result.mtz_reflections.size() > 50, "too few integrated spots");

    const auto oracle = read_oracle_intensities(oracle_path);
    const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_name("I23");
    CHECK_MSG(sg != nullptr, "I23");
    const auto oracle_exp = expand_symmetry(oracle, sg);

    std::vector<double> ratios;
    std::vector<double> abs_z;
    std::size_t n_sigma_ok = 0;
    for (const auto& r : result.mtz_reflections) {
        if (r.sigma <= 0.0f) {
            continue;
        }
        ++n_sigma_ok;
        const auto it = oracle_exp.find({r.h, r.k, r.l});
        if (it == oracle_exp.end()) {
            continue;
        }
        if (it->second.intensity <= 0.0 || r.intensity <= 0.0f ||
            it->second.sigma <= 0.0) {
            continue;
        }
        ratios.push_back(
            static_cast<double>(r.intensity) / it->second.intensity
        );
        abs_z.push_back(std::fabs(
            (static_cast<double>(r.intensity) - it->second.intensity) /
            it->second.sigma
        ));
    }

    CHECK_MSG(ratios.size() >= 50, "too few ASU-matched positive I pairs");
    std::sort(ratios.begin(), ratios.end());
    std::sort(abs_z.begin(), abs_z.end());
    const double med_ratio = ratios[ratios.size() / 2];
    const double med_abs_z = abs_z[abs_z.size() / 2];
    const double frac_ratio_ok = static_cast<double>(std::count_if(
                                     ratios.begin(), ratios.end(),
                                     [](double r) { return r >= 0.1 && r <= 10.0; }
                                 )) /
                                 static_cast<double>(ratios.size());
    const double frac_soft_z =
        static_cast<double>(std::count_if(
            abs_z.begin(), abs_z.end(), [](double z) { return z < 50.0; }
        )) /
        static_cast<double>(abs_z.size());

    std::cout << "insulin_isi: integrated=" << result.mtz_reflections.size()
              << " overlapped=" << result.n_overlapped
              << " matched_pos=" << ratios.size()
              << " med_I_ratio=" << med_ratio << " med_|z|=" << med_abs_z
              << " frac_ratio_[0.1,10]=" << frac_ratio_ok
              << " frac_|z|<50=" << frac_soft_z
              << " sigma_ok=" << n_sigma_ok << '/'
              << result.mtz_reflections.size() << '\n';

    CHECK_MSG(n_sigma_ok == result.mtz_reflections.size(), "non-positive SIGI");
    // Tightened soft gates (E8 PROCESS parity vs F77 oracle).
    CHECK_MSG(med_ratio > 0.3 && med_ratio < 5.0, "median I ratio out of band");
    CHECK_MSG(frac_ratio_ok > 0.55, "too few ratios in [0.1, 10]");
    CHECK_MSG(frac_soft_z > 0.75, "too many |ΔI|/σ_oracle outliers");

    return mosflm_test::finish("insulin_isi");
}
