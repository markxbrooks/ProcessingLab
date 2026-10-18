/**
 * Known-indexed XDET/YDET: C++ PREDICT vs insulin MTZ pixel coordinates.
 *
 * Matches on exact HKL. Applies MOSFLM ``ROTMAT2`` detector normal with
 * refined TILT/TWIST/YSCALE/XCEN/YCEN from the oracle PROCESS log.
 */

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/geometry/LabFrame.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "test_util.hpp"

#include <gemmi/mtz.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
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
    req.detector.distance_mm = 160.0;
    req.detector.beam_x_mm = 0.0;
    req.detector.beam_y_mm = 0.0;
    req.detector.pixel_size_x_mm = 0.0816;
    req.detector.pixel_size_y_mm = 0.0816;
    req.detector.nfast = 2304;
    req.detector.nslow = 2304;
    req.detector.wavelength_A = 0.9790;
    req.detector.yscale = 0.9921;
    req.detector.tilt_mrad = -83.0;
    req.detector.twist_mrad = -140.0;
    req.detector.use_rotmat2 = true;
    req.scan.phi_start_deg = 0.0;
    req.scan.phi_end_deg = 1.0;
    req.settings.dmin_A = 1.45;
    req.settings.ewald_tolerance_Ainv = 0.0025;
    req.settings.mosaic_spread_deg = 0.30;
    req.settings.phi_step_deg = 0.025;
    return req;
}

}  // namespace

int main() {
    const std::filesystem::path mtz =
        "testdata/oracle/known_indexed/insulin_1_001.mtz";
    CHECK_MSG(std::filesystem::exists(mtz), "missing insulin MTZ");

    const auto refs = read_mtz_xdet(mtz);
    CHECK_MSG(refs.size() > 100, "MTZ too small");

    const auto req = insulin_request();
    const auto spots = ReflectionPredictor{}.predict(req).reflections;
    std::map<Hkl, std::pair<double, double>> pred_mm;
    for (const auto& s : spots) {
        const Hkl key{s.h, s.k, s.l};
        if (pred_mm.count(key) == 0U) {
            pred_mm[key] = {s.x_mm, s.y_mm};
        }
    }

    constexpr double kBeamXmm = 94.07;
    constexpr double kBeamYmm = 95.22;
    constexpr double kPixel = 0.0816;
    const double pixel_y = kPixel / req.detector.yscale;

    std::vector<double> residuals;
    residuals.reserve(refs.size());
    for (const auto& ref : refs) {
        const auto it = pred_mm.find({ref.h, ref.k, ref.l});
        if (it == pred_mm.end()) {
            continue;
        }
        double xdet = 0.0;
        double ydet = 0.0;
        mosflm::cpp_mm_to_mosflm_xdet(
            it->second.first,
            it->second.second,
            kBeamXmm,
            kBeamYmm,
            kPixel,
            pixel_y,
            xdet,
            ydet
        );
        residuals.push_back(std::hypot(xdet - ref.xdet, ydet - ref.ydet));
    }

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
              << " inlier_rms_px=" << inlier_rms << '\n';

    // ROTMAT2 + refined TILT/TWIST (was ~40 px flat).
    CHECK_MSG(median < 30.0, "median XDET residual above 30 px");
    CHECK_MSG(inlier_frac > 0.70, "inlier fraction below 0.70");
    CHECK_MSG(inlier_rms < 35.0, "inlier RMS above 35 px");

    return mosflm_test::finish("insulin_xdet");
}
