/**
 * Thaumatin PROFILE 3D PROCESS vs DIALS integrated.mtz (frames 2–10).
 *
 * Soft-skip unless the HDF5 master and DIALS MTZ exist. Gate: Pearson
 * corr(I) > 0.7 on DIALS-strong spots uniquely matched by XDET/ROT
 * (Miller indexing may differ by Friedel / axis permutation).
 */

#include "mosflm/control/ControlSession.hpp"
#include "test_util.hpp"

#include <gemmi/mtz.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

using mosflm::control::ControlSession;

namespace {

struct Obs {
    double i{};
    double ipr{};
    double rot{};
    double xdet{};
    double ydet{};
    double sig{};
};

std::filesystem::path find_master() {
    if (const char* e = std::getenv("MOSFLM_THAUMATIN_H5")) {
        return e;
    }
    if (const char* home = std::getenv("HOME")) {
        return std::filesystem::path(home) / "Downloads" /
               "thaumatin_9_1_master.h5";
    }
    return "thaumatin_9_1_master.h5";
}

std::filesystem::path find_dials_mtz() {
    if (const char* e = std::getenv("MOSFLM_DIALS_MTZ")) {
        return e;
    }
    if (const char* home = std::getenv("HOME")) {
        return std::filesystem::path(home) / "Downloads" /
               "thaumatin_9_1_master" / "2026-10-15-20-29-46-dials" /
               "integrated.mtz";
    }
    return "integrated.mtz";
}

void skip(const std::string& why) {
    std::cout << "thaumatin_profile3d: SKIP (" << why << ")\n";
    std::exit(0);
}

std::vector<Obs> load_dials(const std::filesystem::path& path) {
    gemmi::Mtz m;
    m.read_file(path.string());
    try {
        m.switch_to_original_hkl();
    } catch (const std::exception&) {
    }
    const auto* i = m.column_with_label("I");
    const auto* s = m.column_with_label("SIGI");
    const auto* rot = m.column_with_label("ROT");
    const auto* x = m.column_with_label("XDET");
    const auto* y = m.column_with_label("YDET");
    std::vector<Obs> out;
    if (!i || !s || !rot) {
        return out;
    }
    out.reserve(static_cast<std::size_t>(m.nreflections));
    for (int r = 0; r < m.nreflections; ++r) {
        const std::size_t ii = static_cast<std::size_t>(r);
        Obs o;
        o.i = (*i)[ii];
        o.ipr = o.i;
        o.sig = (*s)[ii];
        o.rot = (*rot)[ii];
        o.xdet = x ? (*x)[ii] : 0.0;
        o.ydet = y ? (*y)[ii] : 0.0;
        out.push_back(o);
    }
    return out;
}

double pearson(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() < 10) {
        return 0.0;
    }
    double mx = 0.0;
    double my = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        mx += x[i];
        my += y[i];
    }
    mx /= static_cast<double>(x.size());
    my /= static_cast<double>(y.size());
    double num = 0.0;
    double dx = 0.0;
    double dy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double ux = x[i] - mx;
        const double uy = y[i] - my;
        num += ux * uy;
        dx += ux * ux;
        dy += uy * uy;
    }
    if (!(dx > 0.0 && dy > 0.0)) {
        return 0.0;
    }
    return num / std::sqrt(dx * dy);
}

}  // namespace

int main() {
#if !defined(MOSFLM_HAS_HDF5)
    skip("built without HDF5");
#endif
    const auto master = find_master();
    const auto dials_path = find_dials_mtz();
    if (!std::filesystem::exists(master)) {
        skip("no " + master.string());
    }
    if (!std::filesystem::exists(dials_path)) {
        skip("no " + dials_path.string());
    }

    ControlSession s;
    s.parse_line("NPROC 1");
    s.parse_line("SYMMETRY 75");
    s.parse_line("WAVELENGTH 0.97949913");
    s.parse_line("DISTANCE 199.78169135");
    s.parse_line("PIXEL 0.075 0.075");
    s.parse_line("BEAM 166.15771234 172.97150533");
    s.parse_line("MOSAIC 0.15");
    s.parse_line("POLARIZATION 0.98");
    s.parse_line("SEPARATION 0.75 0.75");
    s.parse_line("PEAKHALF 5");
    s.parse_line("DMIN 1.15");
    s.parse_line("DMAX 200.0");
    s.parse_line("MINFRAC 0.5");
    s.parse_line("PARTIALSCALE OFF");
    s.parse_line("AUTOSIGMA OFF");
    s.parse_line("POSTREFINE OFF");
    s.parse_line("CELLFIX ON");
    s.parse_line("CELL 57.93557545 57.93557545 150.27023041 90.0 90.0 90.0");
    s.parse_line("PROFILE 3D");
    s.parse_line("NADD 5");
    s.parse_line("MATRIX");
    s.parse_line("-0.0032678601509884 0.0068375384503047 -0.0059789665817106");
    s.parse_line("0.016378382887516 0.0053355240926937 -0.00042364174865803");
    s.parse_line("0.0043584758562717 -0.014923367527643 -0.0028908866964746");
    s.parse_line("PROCESS 2 TO 10");
    const auto result = s.go(master);
    CHECK(result.ran_process);
    CHECK_MSG(
        result.process.mosflm_reflections.size() > 50,
        "too few 3D-integrated reflections"
    );

    std::vector<Obs> pl;
    pl.reserve(result.process.mosflm_reflections.size());
    double rot_max = 0.0;
    for (const auto& row : result.process.mosflm_reflections) {
        Obs o;
        o.i = static_cast<double>(row.intensity);
        o.ipr = static_cast<double>(row.ipr);
        o.sig = static_cast<double>(row.sigma);
        o.rot = static_cast<double>(row.rot);
        o.xdet = static_cast<double>(row.xdet);
        o.ydet = static_cast<double>(row.ydet);
        rot_max = std::max(rot_max, o.rot);
        pl.push_back(o);
    }

    auto dials = load_dials(dials_path);
    CHECK_MSG(!dials.empty(), "failed to read DIALS MTZ");
    std::vector<const Obs*> strong_di;
    for (const auto& d : dials) {
        if (d.rot > rot_max + 0.2 || d.rot < -0.5) {
            continue;
        }
        if (d.sig > 0.0 && d.i / d.sig >= 10.0) {
            strong_di.push_back(&d);
        }
    }

    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> off;
    std::vector<char> used(pl.size(), 0);
    for (const Obs* d : strong_di) {
        const Obs* best = nullptr;
        std::size_t best_i = 0;
        double best_d = 1e99;
        for (std::size_t i = 0; i < pl.size(); ++i) {
            if (used[i]) {
                continue;
            }
            const auto& o = pl[i];
            if (std::fabs(d->rot - o.rot) > 2.0) {
                continue;
            }
            const double dxy = std::hypot(o.xdet - d->xdet, o.ydet - d->ydet);
            if (dxy < best_d) {
                best_d = dxy;
                best = &o;
                best_i = i;
            }
        }
        if (best == nullptr || best_d > 5.0) {
            continue;
        }
        used[best_i] = 1;
        a.push_back(best->i);
        b.push_back(d->i);
        off.push_back(best_d);
    }

    double med_off = 0.0;
    if (!off.empty()) {
        std::nth_element(
            off.begin(), off.begin() + static_cast<std::ptrdiff_t>(off.size() / 2),
            off.end()
        );
        med_off = off[off.size() / 2];
    }
    const double cc = pearson(a, b);
    std::cout << "thaumatin_profile3d: nref="
              << result.process.mosflm_reflections.size()
              << " matched_strong=" << a.size() << " cc(I)=" << cc
              << " med|dxy|=" << med_off << " n_dials_strong="
              << strong_di.size() << '\n';
    CHECK_MSG(a.size() >= 40, "too few strong matched spots vs DIALS");
    CHECK_MSG(cc > 0.7, "thaumatin 3D corr(I) vs DIALS <= 0.7");

    return mosflm_test::finish("thaumatin_profile3d");
}
