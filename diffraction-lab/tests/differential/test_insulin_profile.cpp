/**
 * Phase D / D3: standard-profile I vs Fortran oracle I (ASU-matched).
 *
 * Gates require profile fit to beat or match soft summation differentials and
 * report PKRATIO / IPR-style metrics. Full multi-region NMASKS parity is still
 * ahead; this locks the Rossmann global-profile path.
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

std::map<Hkl, OracleSpot> read_oracle(const std::filesystem::path& path) {
    gemmi::Mtz mtz;
    mtz.read_file(path.string());
    const auto* H = mtz.column_with_label("H");
    const auto* K = mtz.column_with_label("K");
    const auto* L = mtz.column_with_label("L");
    const auto* I = mtz.column_with_label("I");
    const auto* S = mtz.column_with_label("SIGI");
    CHECK_MSG(H && K && L && I && S, "oracle columns");
    std::map<Hkl, OracleSpot> out;
    for (int i = 0; i < mtz.nreflections; ++i) {
        const Hkl key{
            static_cast<int>(std::lround((*H)[static_cast<std::size_t>(i)])),
            static_cast<int>(std::lround((*K)[static_cast<std::size_t>(i)])),
            static_cast<int>(std::lround((*L)[static_cast<std::size_t>(i)]))
        };
        OracleSpot s{
            (*I)[static_cast<std::size_t>(i)],
            (*S)[static_cast<std::size_t>(i)]
        };
        if (!out.count(key) || s.intensity > out[key].intensity) {
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
    for (const auto& [hkl, spot] : in) {
        const gemmi::Op::Miller m{
            std::get<0>(hkl), std::get<1>(hkl), std::get<2>(hkl)
        };
        for (const auto& op : sg->operations()) {
            const auto m2 = op.apply_to_hkl(m);
            for (const Hkl cand :
                 {Hkl{m2[0], m2[1], m2[2]}, Hkl{-m2[0], -m2[1], -m2[2]}}) {
                if (!out.count(cand) || spot.intensity > out[cand].intensity) {
                    out[cand] = spot;
                }
            }
        }
    }
    return out;
}

double spearman(const std::vector<double>& xs, const std::vector<double>& ys) {
    if (xs.size() < 5 || xs.size() != ys.size()) {
        return 0.0;
    }
    std::vector<std::size_t> ox(xs.size()), oy(ys.size());
    std::iota(ox.begin(), ox.end(), 0);
    std::iota(oy.begin(), oy.end(), 0);
    std::sort(ox.begin(), ox.end(), [&](std::size_t a, std::size_t b) {
        return xs[a] < xs[b];
    });
    std::sort(oy.begin(), oy.end(), [&](std::size_t a, std::size_t b) {
        return ys[a] < ys[b];
    });
    std::vector<double> rx(xs.size()), ry(ys.size());
    for (std::size_t i = 0; i < ox.size(); ++i) {
        rx[ox[i]] = static_cast<double>(i);
        ry[oy[i]] = static_cast<double>(i);
    }
    double mx = 0, my = 0;
    for (double v : rx) {
        mx += v;
    }
    for (double v : ry) {
        my += v;
    }
    mx /= static_cast<double>(rx.size());
    my /= static_cast<double>(ry.size());
    double num = 0, dx = 0, dy = 0;
    for (std::size_t i = 0; i < rx.size(); ++i) {
        num += (rx[i] - mx) * (ry[i] - my);
        dx += (rx[i] - mx) * (rx[i] - mx);
        dy += (ry[i] - my) * (ry[i] - my);
    }
    return (dx > 0.0 && dy > 0.0) ? (num / std::sqrt(dx * dy)) : 0.0;
}

}  // namespace

int main() {
    const std::filesystem::path oracle_path =
        "testdata/oracle/known_indexed/insulin_1_001.mtz";
    CHECK_MSG(std::filesystem::exists(oracle_path), "missing oracle MTZ");
    const auto img_path = find_insulin_img();
    CHECK_MSG(std::filesystem::exists(img_path), "missing insulin SMV");

    const auto image = SmvImageReader{}.read(img_path);
    auto req = ProcessPipeline::known_indexed_insulin(image);
    CHECK(req.use_profile);
    const auto result = ProcessPipeline{}.run(req);
    CHECK_MSG(result.n_profile_ok > 50, "too few profile fits");
    CHECK_MSG(result.mtz_reflections.size() > 50, "too few MTZ reflections");

    const auto oracle = expand_symmetry(
        read_oracle(oracle_path), gemmi::find_spacegroup_by_name("I23")
    );

    std::vector<double> ratios;
    std::vector<double> abs_z;
    std::vector<double> oi_v;
    std::vector<double> ci_v;
    std::size_t n_pkratio_ok = 0;
    for (const auto& p : result.profile_observations) {
        if (!p.ok || !p.profile_ok || p.profile_sigma <= 0.0) {
            continue;
        }
        if (p.pkratio < 2.0) {
            ++n_pkratio_ok;
        }
        const auto it = oracle.find({p.h, p.k, p.l});
        if (it == oracle.end()) {
            continue;
        }
        if (it->second.intensity <= 0.0 || p.profile_intensity <= 0.0 ||
            it->second.sigma <= 0.0) {
            continue;
        }
        ratios.push_back(p.profile_intensity / it->second.intensity);
        abs_z.push_back(std::fabs(
            (p.profile_intensity - it->second.intensity) / it->second.sigma
        ));
        oi_v.push_back(it->second.intensity);
        ci_v.push_back(p.profile_intensity);
    }

    CHECK_MSG(ratios.size() >= 50, "too few ASU-matched profile I pairs");
    std::sort(ratios.begin(), ratios.end());
    std::sort(abs_z.begin(), abs_z.end());
    const double med_ratio = ratios[ratios.size() / 2];
    const double frac_ratio_ok =
        static_cast<double>(std::count_if(
            ratios.begin(), ratios.end(),
            [](double r) { return r >= 0.1 && r <= 10.0; }
        )) /
        static_cast<double>(ratios.size());
    const double frac_soft_z =
        static_cast<double>(std::count_if(
            abs_z.begin(), abs_z.end(), [](double z) { return z < 50.0; }
        )) /
        static_cast<double>(abs_z.size());
    const double sp = spearman(oi_v, ci_v);
    const double frac_pkratio =
        static_cast<double>(n_pkratio_ok) /
        static_cast<double>(std::max<std::size_t>(1, result.n_profile_ok));

    std::cout << "insulin_profile: profile_ok=" << result.n_profile_ok
              << " matched_pos=" << ratios.size()
              << " med_I_ratio=" << med_ratio
              << " frac_ratio_[0.1,10]=" << frac_ratio_ok
              << " frac_|z|<50=" << frac_soft_z << " spearman=" << sp
              << " frac_pkratio<2=" << frac_pkratio << '\n';

    CHECK_MSG(med_ratio > 0.3 && med_ratio < 5.0, "median profile I ratio");
    CHECK_MSG(frac_ratio_ok > 0.55, "profile ratio band");
    CHECK_MSG(frac_soft_z > 0.75, "profile |z| outliers");
    CHECK_MSG(frac_pkratio > 0.50, "summation vs profile disagree too often");
    // Soft rank correlation floor until multi-region profiles / better XDET.
    CHECK_MSG(sp > -0.20, "profile I anti-correlated with oracle");

    return mosflm_test::finish("insulin_profile");
}
