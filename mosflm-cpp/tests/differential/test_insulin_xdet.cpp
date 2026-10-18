/**
 * Known-indexed XDET/YDET: C++ PREDICT vs insulin MTZ pixel coordinates.
 *
 * Matches on exact HKL (ASU choice differs for many spots). Soft residual
 * gates while TILT/TWIST/YSCALE are not modelled in DetectorGeometry.
 */

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/geometry/LabFrame.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "test_util.hpp"

#include <gemmi/mtz.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

using mosflm::OrientationMatrix;
using mosflm::crystal::CrystalMetadata;
using mosflm::diffraction::PredictionRequest;
using mosflm::diffraction::ReflectionPredictor;

namespace {

using Hkl = std::tuple<int, int, int>;

struct MtzSpot {
    int h{};
    int k{};
    int l{};
    double xdet{};
    double ydet{};
};

std::vector<MtzSpot> read_mtz_xdet(const std::filesystem::path& path) {
    gemmi::Mtz mtz;
    mtz.read_file(path.string());
    const gemmi::Mtz::Column* col_h = mtz.column_with_label("H");
    const gemmi::Mtz::Column* col_k = mtz.column_with_label("K");
    const gemmi::Mtz::Column* col_l = mtz.column_with_label("L");
    const gemmi::Mtz::Column* col_x = mtz.column_with_label("XDET");
    const gemmi::Mtz::Column* col_y = mtz.column_with_label("YDET");
    if (!col_h || !col_k || !col_l || !col_x || !col_y) {
        throw std::runtime_error("MTZ missing H/K/L/XDET/YDET columns");
    }
    std::vector<MtzSpot> out;
    out.reserve(static_cast<std::size_t>(mtz.nreflections));
    for (int i = 0; i < mtz.nreflections; ++i) {
        MtzSpot s;
        s.h = static_cast<int>(std::lround((*col_h)[static_cast<std::size_t>(i)]));
        s.k = static_cast<int>(std::lround((*col_k)[static_cast<std::size_t>(i)]));
        s.l = static_cast<int>(std::lround((*col_l)[static_cast<std::size_t>(i)]));
        s.xdet = (*col_x)[static_cast<std::size_t>(i)];
        s.ydet = (*col_y)[static_cast<std::size_t>(i)];
        out.push_back(s);
    }
    return out;
}

PredictionRequest insulin_request() {
    Eigen::Matrix3d A;
    A << -0.012237, -0.000358, -0.002626, 0.002191, 0.005616, -0.010974,
        0.001492, -0.011185, -0.005426;

    PredictionRequest req;
    req.crystal = CrystalMetadata::from_parameters(
        78.190, 78.190, 78.190, 90.0, 90.0, 90.0, "I 2 3"
    );
    req.orientation = OrientationMatrix::from_mosflm_orientation_matrix(A);
    if (std::getenv("LAMSCALE") != nullptr) {
        req.orientation = OrientationMatrix(req.orientation.ub() / 0.9790);
    }
    req.detector.distance_mm = 160.0;
    req.detector.beam_x_mm = 0.0;
    req.detector.beam_y_mm = 0.0;
    req.detector.pixel_size_x_mm = 0.0816;
    req.detector.pixel_size_y_mm = 0.0816;
    req.detector.nfast = 2304;
    req.detector.nslow = 2304;
    req.detector.wavelength_A = 0.9790;
    req.scan.phi_start_deg = 0.0;
    req.scan.phi_end_deg = 1.0;
    req.settings.dmin_A = 1.45;
    req.settings.ewald_tolerance_Ainv = 0.0025;
    req.settings.phi_step_deg = 0.025;
    return req;
}

Eigen::Matrix3d rot_x_deg(double phi_deg) {
    const double r = phi_deg * (3.14159265358979323846 / 180.0);
    const double c = std::cos(r);
    const double s = std::sin(r);
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
    R(1, 1) = c;
    R(1, 2) = -s;
    R(2, 1) = s;
    R(2, 2) = c;
    return R;
}

/// Exact-HKL matched MTZ spots with the Ewald-satisfying diffracted ray s1.
struct Target {
    MtzSpot mtz;
    Eigen::Vector3d s1;
};

std::vector<Target> build_targets(
    const std::vector<MtzSpot>& refs,
    const PredictionRequest& req
) {
    const double inv_lambda = 1.0 / req.detector.wavelength_A;
    const Eigen::Vector3d s0(0.0, 0.0, inv_lambda);
    const Eigen::Matrix3d UB = req.orientation.ub();
    const double tol = req.settings.ewald_tolerance_Ainv;

    std::set<Hkl> predicted_hkl;
    for (const auto& s : ReflectionPredictor{}.predict(req).reflections) {
        predicted_hkl.insert(Hkl{s.h, s.k, s.l});
    }

    std::vector<Target> targets;
    for (const auto& m : refs) {
        if (predicted_hkl.find(Hkl{m.h, m.k, m.l}) == predicted_hkl.end()) {
            continue;
        }
        const Eigen::Vector3d q0 = UB * Eigen::Vector3d(m.h, m.k, m.l);
        Eigen::Vector3d s1;
        double best = 1e18;
        for (double phi = 0.0; phi <= 1.0 + 1e-9; phi += 0.01) {
            const Eigen::Vector3d v = s0 + rot_x_deg(phi) * q0;
            const double e = std::fabs(v.norm() - inv_lambda);
            if (e < best) {
                best = e;
                s1 = v;
            }
        }
        if (best <= tol) {
            targets.push_back({m, s1});
        }
    }
    return targets;
}

double median_of(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return v.empty() ? 0.0 : v[v.size() / 2];
}

void project(
    const mosflm::DetectorGeometry& det,
    const Eigen::Vector3d& s1,
    double& xdet,
    double& ydet
) {
    double x_mm = 0.0;
    double y_mm = 0.0;
    if (!det.lab_to_detector_mm(s1, x_mm, y_mm)) {
        xdet = std::nan("");
        ydet = std::nan("");
        return;
    }
    mosflm::cpp_mm_to_mosflm_xdet(
        x_mm, y_mm, 94.0, 94.7, 0.0816, 0.0816, xdet, ydet
    );
}

std::vector<double> residuals_for(
    const PredictionRequest& req,
    const std::vector<MtzSpot>& refs,
    double bx,
    double by
) {
    const auto spots = ReflectionPredictor{}.predict(req).reflections;
    std::map<Hkl, std::pair<double, double>> pred;
    for (const auto& s : spots) {
        pred[Hkl{s.h, s.k, s.l}] = {s.x_mm, s.y_mm};
    }
    std::vector<double> residuals;
    for (const auto& m : refs) {
        const auto it = pred.find(Hkl{m.h, m.k, m.l});
        if (it == pred.end()) {
            continue;
        }
        double xdet = 0.0;
        double ydet = 0.0;
        mosflm::cpp_mm_to_mosflm_xdet(
            it->second.first, it->second.second, bx, by, 0.0816, 0.0816, xdet,
            ydet
        );
        residuals.push_back(std::hypot(xdet - m.xdet, ydet - m.ydet));
    }
    return residuals;
}

void dump(const std::vector<MtzSpot>& refs) {
    const auto req = insulin_request();
    const auto targets = build_targets(refs, req);
    std::printf("T count=%zu\n", targets.size());
    for (const auto& t : targets) {
        double x_mm = 0.0;
        double y_mm = 0.0;
        if (!req.detector.lab_to_detector_mm(t.s1, x_mm, y_mm)) {
            continue;
        }
        std::printf(
            "T %d %d %d %.6f %.6f %.4f %.4f\n", t.mtz.h, t.mtz.k, t.mtz.l, x_mm,
            y_mm, t.mtz.xdet, t.mtz.ydet
        );
    }
}

void sweep(const std::vector<MtzSpot>& refs) {
    const auto req = insulin_request();
    const auto targets = build_targets(refs, req);
    std::printf("SWEEP targets=%zu\n", targets.size());

    auto eval = [&](double tilt, double twist, double ys) {
        auto det = req.detector;
        det.tilt_mrad = tilt;
        det.twist_mrad = twist;
        det.yscale = ys;
        std::vector<double> r;
        r.reserve(targets.size());
        for (const auto& t : targets) {
            double xdet = 0.0;
            double ydet = 0.0;
            project(det, t.s1, xdet, ydet);
            if (std::isnan(xdet)) {
                continue;
            }
            r.push_back(std::hypot(xdet - t.mtz.xdet, ydet - t.mtz.ydet));
        }
        return median_of(r);
    };

    std::printf(
        "SWEEP header median=%.3f\n", eval(-10.996, -18.675, 0.992076)
    );

    double best_med = 1e18;
    double bt = 0.0;
    double bw = 0.0;
    double by = 1.0;
    for (double tilt = -60.0; tilt <= 60.0; tilt += 1.5) {
        for (double twist = -60.0; twist <= 60.0; twist += 1.5) {
            for (double ys = 0.985; ys <= 1.0001; ys += 0.001) {
                const double med = eval(tilt, twist, ys);
                if (med < best_med) {
                    best_med = med;
                    bt = tilt;
                    bw = twist;
                    by = ys;
                }
            }
        }
    }
    std::printf(
        "SWEEP best median=%.3f px at tilt_mrad=%.1f twist_mrad=%.1f "
        "yscale=%.4f\n",
        best_med, bt, bw, by
    );
}

}  // namespace

int main() {
    const std::filesystem::path mtz =
        "testdata/oracle/known_indexed/insulin_1_001.mtz";
    CHECK_MSG(std::filesystem::exists(mtz), "missing insulin MTZ");

    const auto refs = read_mtz_xdet(mtz);
    CHECK_MSG(refs.size() > 100, "MTZ too small");

    if (std::getenv("SWEEP") != nullptr) {
        sweep(refs);
        return 0;
    }
    if (std::getenv("DUMP") != nullptr) {
        dump(refs);
        return 0;
    }

    auto req = insulin_request();
    if (const char* t = std::getenv("TILT")) {
        req.detector.tilt_mrad = std::atof(t);
    }
    if (const char* w = std::getenv("TWIST")) {
        req.detector.twist_mrad = std::atof(w);
    }
    if (const char* y = std::getenv("YSCALE")) {
        req.detector.yscale = std::atof(y);
    }
    double bx = 94.0;
    double by = 94.7;
    if (const char* b = std::getenv("BX")) {
        bx = std::atof(b);
    }
    if (const char* b = std::getenv("BY")) {
        by = std::atof(b);
    }

    auto residuals = residuals_for(req, refs, bx, by);

    CHECK_MSG(residuals.size() > 50, "too few exact-HKL positional matches");
    std::sort(residuals.begin(), residuals.end());
    const double median = residuals[residuals.size() / 2];
    const double p90 = residuals[(9 * residuals.size()) / 10];

    std::size_t n_inlier = 0;
    double inlier_sum2 = 0.0;
    for (double r : residuals) {
        if (r < 100.0) {
            ++n_inlier;
            inlier_sum2 += r * r;
        }
    }
    const double inlier_frac =
        static_cast<double>(n_inlier) / static_cast<double>(residuals.size());
    const double inlier_rms =
        n_inlier > 0 ? std::sqrt(inlier_sum2 / static_cast<double>(n_inlier))
                     : 0.0;

    std::cout << "insulin_xdet: exact_matches=" << residuals.size()
              << " median_px=" << median << " p90_px=" << p90
              << " inlier<100frac=" << inlier_frac
              << " inlier_rms_px=" << inlier_rms
              << " tilt_mrad=" << req.detector.tilt_mrad
              << " twist_mrad=" << req.detector.twist_mrad
              << " yscale=" << req.detector.yscale << '\n';

    return 0;
}
