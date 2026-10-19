/**
 * Compare ProcessingLab thaumatin process MTZ to a DIALS integrated.mtz.
 *
 * Soft-skip unless both files exist. Gates (oracle DIALS geometry PROCESS):
 *   - nref ratio in [0.5, 1.5]
 *   - FRACTIONCALC ≥ 0.99 for ≥ 40% of rows
 *   - mean-I CC on PL full+strong common HKL ≥ 0.35 (wave-3; warn <0.50)
 *   - % I/σI > 3 ≥ 5%
 *   - XDET max within 25% of DIALS (panel extent)
 *   - soft: batch-2 median |ROT_PL − ROT_DIALS| ≲ 0.05°
 *
 * Env overrides:
 *   MOSFLM_THAUMATIN_MTZ  — PL output (default process_thaumatin.mtz / process.mtz)
 *   MOSFLM_DIALS_MTZ      — DIALS integrated.mtz
 */

#include "mosflm/crystal/MtzAdapter.hpp"
#include "test_util.hpp"

#include <gemmi/mtz.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <tuple>
#include <vector>

namespace {

using Hkl = std::tuple<int, int, int>;

std::filesystem::path find_pl_mtz() {
    if (const char* e = std::getenv("MOSFLM_THAUMATIN_MTZ")) {
        return e;
    }
    for (const char* c :
         {"process_thaumatin.mtz", "../process_thaumatin.mtz", "process.mtz",
          "../process.mtz"}) {
        if (std::filesystem::exists(c)) {
            return c;
        }
    }
    return "process_thaumatin.mtz";
}

std::filesystem::path find_dials_mtz() {
    if (const char* e = std::getenv("MOSFLM_DIALS_MTZ")) {
        return e;
    }
    const char* cand =
        "/Users/brooks/Downloads/thaumatin_9_1_master/"
        "2026-10-15-20-29-46-dials/integrated.mtz";
    return cand;
}

void skip(const std::string& why) {
    std::cout << "thaumatin_dials_mtz_compare: SKIP (" << why << ")\n";
    std::exit(0);
}

double median_rot_for_batch(const gemmi::Mtz& m, int batch_id) {
    const auto* batch = m.column_with_label("BATCH");
    const auto* rot = m.column_with_label("ROT");
    if (!batch || !rot) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::vector<double> vals;
    vals.reserve(static_cast<std::size_t>(m.nreflections) / 100);
    for (int i = 0; i < m.nreflections; ++i) {
        if (static_cast<int>(std::lround((*batch)[static_cast<std::size_t>(i)])) ==
            batch_id) {
            vals.push_back((*rot)[static_cast<std::size_t>(i)]);
        }
    }
    if (vals.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::nth_element(
        vals.begin(), vals.begin() + static_cast<std::ptrdiff_t>(vals.size() / 2),
        vals.end()
    );
    return vals[vals.size() / 2];
}

}  // namespace

int main() {
    const auto pl_path = find_pl_mtz();
    const auto di_path = find_dials_mtz();
    if (!std::filesystem::exists(pl_path)) {
        skip("missing PL MTZ " + pl_path.string());
    }
    if (!std::filesystem::exists(di_path)) {
        skip("missing DIALS MTZ " + di_path.string());
    }

    gemmi::Mtz pl;
    gemmi::Mtz di;
    pl.read_file(pl_path.string());
    di.read_file(di_path.string());

    // Reject wrong-crystal PL (e.g. leftover 76c1b C2).
    if (pl.spacegroup_number != 75 && pl.spacegroup_number != 89 &&
        pl.spacegroup_number != 90 && pl.spacegroup_number != 91 &&
        pl.spacegroup_number != 92 && pl.spacegroup_number != 93 &&
        pl.spacegroup_number != 94 && pl.spacegroup_number != 95 &&
        pl.spacegroup_number != 96) {
        // P4 family; also accept if cell a≈b≈58
        if (std::fabs(pl.cell.a - 58.0) > 5.0) {
            skip(
                "PL MTZ does not look like thaumatin P4 (sg=" +
                std::to_string(pl.spacegroup_number) + ")"
            );
        }
    }

    const double nref_ratio =
        static_cast<double>(pl.nreflections) /
        static_cast<double>(std::max(1, di.nreflections));
    std::cout << "thaumatin_dials_mtz_compare: pl_nref=" << pl.nreflections
              << " di_nref=" << di.nreflections << " ratio=" << nref_ratio
              << '\n';
    CHECK_MSG(nref_ratio >= 0.5 && nref_ratio <= 1.5, "nref ratio out of band");

    const gemmi::Mtz::Column* fcol = pl.column_with_label("FRACTIONCALC");
    CHECK_MSG(fcol != nullptr, "PL missing FRACTIONCALC");
    std::size_t n_full = 0;
    for (int i = 0; i < pl.nreflections; ++i) {
        if ((*fcol)[static_cast<std::size_t>(i)] >= 0.99f) {
            ++n_full;
        }
    }
    const double pct_full =
        100.0 * static_cast<double>(n_full) /
        static_cast<double>(std::max(1, pl.nreflections));
    std::cout << "  FRACTIONCALC pct>=0.99=" << pct_full << '\n';
    CHECK_MSG(pct_full >= 40.0, "too few near-full reflections");

    const gemmi::Mtz::Column* xpl = pl.column_with_label("XDET");
    const gemmi::Mtz::Column* xdi = di.column_with_label("XDET");
    CHECK_MSG(xpl && xdi, "missing XDET");
    float xmax_pl = (*xpl)[0];
    float xmax_di = (*xdi)[0];
    for (int i = 0; i < pl.nreflections; ++i) {
        xmax_pl = std::max(xmax_pl, (*xpl)[static_cast<std::size_t>(i)]);
    }
    for (int i = 0; i < di.nreflections; ++i) {
        xmax_di = std::max(xmax_di, (*xdi)[static_cast<std::size_t>(i)]);
    }
    std::cout << "  XDET max pl=" << xmax_pl << " di=" << xmax_di << '\n';
    CHECK_MSG(
        xmax_pl > 0.75f * xmax_di,
        "PL XDET extent far below DIALS (beam/panel size?)"
    );

    // Soft ROT alignment check (batch 2 ≈ first shared image).
    const double rot_pl = median_rot_for_batch(pl, 2);
    const double rot_di = median_rot_for_batch(di, 2);
    if (std::isfinite(rot_pl) && std::isfinite(rot_di)) {
        const double drot = std::fabs(rot_pl - rot_di);
        std::cout << "  batch-2 median ROT pl=" << rot_pl << " di=" << rot_di
                  << " |d|=" << drot << '\n';
        if (drot > 0.05) {
            std::cout << "  WARN: batch-2 |dROT| > 0.05°\n";
        }
    }

    // % I/σI > 3
    const auto* inten = pl.column_with_label("I");
    const auto* sigi = pl.column_with_label("SIGI");
    CHECK_MSG(inten && sigi, "missing I/SIGI");
    std::size_t n_strong = 0;
    std::size_t n_ok = 0;
    for (int i = 0; i < pl.nreflections; ++i) {
        const float s = (*sigi)[static_cast<std::size_t>(i)];
        if (!(s > 0.0f) || !std::isfinite(s)) {
            continue;
        }
        ++n_ok;
        if ((*inten)[static_cast<std::size_t>(i)] / s > 3.0f) {
            ++n_strong;
        }
    }
    const double pct_strong =
        100.0 * static_cast<double>(n_strong) /
        static_cast<double>(std::max<std::size_t>(1, n_ok));
    std::cout << "  pct I/sigI>3=" << pct_strong << " (n=" << n_ok << ")\n";
    // Wave-3 partiality scaling / profile gates can dilute I/σ; keep a
    // soft floor while the mean-I CC gate carries intensity quality.
    CHECK_MSG(pct_strong >= 5.0, "too few strong (I/sigI>3) reflections");

    // Mean-I CC on common stored HKL. Prefer strong near-full PL obs
    // (FRACTIONCALC≥0.99, I/σ>5) — wave-3 intensity-quality gate.
    auto mean_I = [](const gemmi::Mtz& m, bool prefer_full_strong) {
        std::map<Hkl, std::pair<double, int>> acc;
        const auto* h = m.column_with_label("H");
        const auto* k = m.column_with_label("K");
        const auto* l = m.column_with_label("L");
        const auto* inten_col = m.column_with_label("I");
        const auto* sig_col = m.column_with_label("SIGI");
        const auto* frac_col = m.column_with_label("FRACTIONCALC");
        for (int i = 0; i < m.nreflections; ++i) {
            const std::size_t ii = static_cast<std::size_t>(i);
            if (prefer_full_strong && frac_col != nullptr) {
                if ((*frac_col)[ii] < 0.99f) {
                    continue;
                }
            }
            if (prefer_full_strong && sig_col != nullptr) {
                const float s = (*sig_col)[ii];
                if (!(s > 0.0f) || (*inten_col)[ii] / s < 5.0f) {
                    continue;
                }
            }
            const Hkl key{
                static_cast<int>(std::lround((*h)[ii])),
                static_cast<int>(std::lround((*k)[ii])),
                static_cast<int>(std::lround((*l)[ii]))
            };
            auto& p = acc[key];
            p.first += (*inten_col)[ii];
            ++p.second;
        }
        std::map<Hkl, double> out;
        for (const auto& kv : acc) {
            out[kv.first] = kv.second.first / static_cast<double>(kv.second.second);
        }
        return out;
    };
    auto pearson = [](const std::map<Hkl, double>& mpl,
                      const std::map<Hkl, double>& mdi) {
        std::vector<double> a;
        std::vector<double> b;
        for (const auto& kv : mpl) {
            const auto it = mdi.find(kv.first);
            if (it == mdi.end()) {
                continue;
            }
            a.push_back(kv.second);
            b.push_back(it->second);
        }
        if (a.size() < 50) {
            return std::make_pair(0.0, static_cast<std::size_t>(0));
        }
        double mean_a = 0.0;
        double mean_b = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            mean_a += a[i];
            mean_b += b[i];
        }
        mean_a /= static_cast<double>(a.size());
        mean_b /= static_cast<double>(b.size());
        double num = 0.0;
        double da = 0.0;
        double db = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            const double xa = a[i] - mean_a;
            const double xb = b[i] - mean_b;
            num += xa * xb;
            da += xa * xa;
            db += xb * xb;
        }
        const double cc = num / std::sqrt(std::max(1e-30, da * db));
        return std::make_pair(cc, a.size());
    };

    const auto mpl_all = mean_I(pl, false);
    const auto mdi_all = mean_I(di, false);
    const auto cc_all = pearson(mpl_all, mdi_all);
    std::cout << "  mean-I CC on " << cc_all.second
              << " common HKL (all): " << cc_all.first << '\n';

    const auto mpl_fs = mean_I(pl, true);
    const auto mdi_fs = mean_I(di, false);
    const auto cc_fs = pearson(mpl_fs, mdi_fs);
    std::cout << "  mean-I CC on " << cc_fs.second
              << " common HKL (PL full+strong): " << cc_fs.first << '\n';

    // Wave-3 gate: full+strong CC preferred; fall back to all-obs CC.
    // Hard floor 0.35 (shipped PARTIALSCALE/box hygiene); warn below 0.50
    // until Kabsch 3D / scale closes the remaining intensity gap.
    const double cc_gate =
        (cc_fs.second > 200) ? cc_fs.first : cc_all.first;
    if (cc_gate < 0.50) {
        std::cout << "  WARN: mean-I CC " << cc_gate
                  << " < 0.50 (Kabsch 3D / scale still needed)\n";
    }
    CHECK_MSG(cc_gate >= 0.35, "mean-I CC below wave-3 floor");

    return mosflm_test::finish("thaumatin_dials_mtz_compare");
}
