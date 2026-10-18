#include "mosflm/diffraction/ReflectionPredictor.hpp"

#include "mosflm/geometry/ReciprocalLattice.hpp"
#include "mosflm/geometry/UnitCell.hpp"

#include <gemmi/symmetry.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mosflm::diffraction {
namespace {

Eigen::Matrix3d rotation_x_deg(double phi_deg) {
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

bool on_detector(const DetectorGeometry& det, double x_mm, double y_mm) {
    double xf = 0.0;
    double yf = 0.0;
    det.mm_to_pixel(x_mm, y_mm, xf, yf);
    return std::fabs(xf) < 0.5 * static_cast<double>(det.nfast) &&
           std::fabs(yf) < 0.5 * static_cast<double>(det.nslow);
}

struct EvalPhi {
    bool ok{false};
    double x_mm{0.0};
    double y_mm{0.0};
};

EvalPhi evaluate_phi(
    const Eigen::Vector3d& q0,
    double phi_deg,
    const Eigen::Vector3d& s0,
    double inv_lambda,
    double ewald_tol,
    const DetectorGeometry& det
) {
    EvalPhi out;
    const Eigen::Vector3d q = rotation_x_deg(phi_deg) * q0;
    const Eigen::Vector3d s1 = s0 + q;
    if (std::fabs(s1.norm() - inv_lambda) > ewald_tol) {
        return out;
    }
    if (!det.lab_to_detector_mm(s1, out.x_mm, out.y_mm)) {
        return out;
    }
    if (!on_detector(det, out.x_mm, out.y_mm)) {
        return out;
    }
    out.ok = true;
    return out;
}

PredictionRequest base_analytical_fixture() {
    PredictionRequest req;
    const UnitCell cell{58.0, 58.0, 150.0, 90.0, 90.0, 90.0};
    req.crystal = crystal::CrystalMetadata::from_mosflm_cell(cell, "P 1");
    req.orientation = OrientationMatrix::from_cell_and_rotation(cell);
    req.detector.wavelength_A = 1.0;
    req.detector.distance_mm = 250.0;
    req.detector.beam_x_mm = 0.0;
    req.detector.beam_y_mm = 0.0;
    req.detector.pixel_size_x_mm = 0.172;
    req.detector.pixel_size_y_mm = 0.172;
    req.detector.nfast = 1475;
    req.detector.nslow = 1679;
    req.settings.ewald_tolerance_Ainv = 0.0025;
    req.settings.phi_step_deg = 0.05;
    return req;
}

}  // namespace

PredictionRequest ReflectionPredictor::analytical_still_fixture() {
    auto req = base_analytical_fixture();
    req.scan.phi_start_deg = 0.0;
    req.scan.phi_end_deg = 0.0;
    req.settings.dmin_A = 3.5;
    return req;
}

PredictionRequest ReflectionPredictor::analytical_oscillation_fixture() {
    auto req = base_analytical_fixture();
    req.scan.phi_start_deg = 0.0;
    req.scan.phi_end_deg = 0.5;
    req.settings.dmin_A = 5.0;
    req.settings.phi_step_deg = 0.05;
    return req;
}

PredictionResult ReflectionPredictor::predict(
    const PredictionRequest& request
) const {
    const double lambda = request.detector.wavelength_A;
    if (!(lambda > 0.0) || !(request.settings.dmin_A > 0.0)) {
        throw std::invalid_argument(
            "ReflectionPredictor::predict: invalid wavelength or dmin"
        );
    }
    if (request.scan.phi_end_deg < request.scan.phi_start_deg) {
        throw std::invalid_argument(
            "ReflectionPredictor::predict: phi_end < phi_start"
        );
    }

    const UnitCell cell = request.crystal.to_mosflm_cell();
    const ReciprocalLattice recip = ReciprocalLattice::from_direct(cell);
    const double inv_lambda = 1.0 / lambda;
    const Eigen::Vector3d s0(0.0, 0.0, inv_lambda);
    const Eigen::Matrix3d& UB = request.orientation.ub();
    const bool still = std::fabs(
                           request.scan.phi_end_deg - request.scan.phi_start_deg
                       ) < 1e-12;

    // Mosaic spread widens the Ewald shell (Å⁻¹ ≈ mosaic_rad / λ).
    const double mosaic_rad =
        request.settings.mosaic_spread_deg * (3.14159265358979323846 / 180.0);
    const double mosaic_tol = mosaic_rad * inv_lambda;
    const double ewald_tol = std::max(
        request.settings.ewald_tolerance_Ainv, mosaic_tol
    );

    gemmi::GroupOps gops;
    const bool filter_abs =
        request.settings.apply_systematic_absences &&
        request.crystal.spacegroup != nullptr;
    if (filter_abs) {
        gops = request.crystal.spacegroup->operations();
    }

    const int max_h =
        static_cast<int>(std::ceil(cell.a / request.settings.dmin_A)) + 1;
    const int max_k =
        static_cast<int>(std::ceil(cell.b / request.settings.dmin_A)) + 1;
    const int max_l =
        static_cast<int>(std::ceil(cell.c / request.settings.dmin_A)) + 1;

    // Build φ samples (inclusive endpoints).
    std::vector<double> phis;
    if (still) {
        phis.push_back(request.scan.phi_start_deg);
    } else {
        const double step = std::max(request.settings.phi_step_deg, 1e-6);
        for (double phi = request.scan.phi_start_deg;
             phi <= request.scan.phi_end_deg + 1e-9; phi += step) {
            phis.push_back(std::min(phi, request.scan.phi_end_deg));
        }
        if (phis.empty() ||
            std::fabs(phis.back() - request.scan.phi_end_deg) > 1e-9) {
            phis.push_back(request.scan.phi_end_deg);
        }
    }

    std::vector<PredictedReflection> spots;
    for (int h = -max_h; h <= max_h; ++h) {
        for (int k = -max_k; k <= max_k; ++k) {
            for (int l = -max_l; l <= max_l; ++l) {
                if (h == 0 && k == 0 && l == 0) {
                    continue;
                }
                if (filter_abs &&
                    gops.is_systematically_absent({h, k, l})) {
                    continue;
                }
                if (recip.d_spacing(h, k, l) < request.settings.dmin_A) {
                    continue;
                }
                const Eigen::Vector3d q0 = UB * Eigen::Vector3d(h, k, l);

                double phi_first = 0.0;
                double phi_last = 0.0;
                double x_mid = 0.0;
                double y_mid = 0.0;
                int n_on = 0;
                bool first_on = false;
                bool last_on = false;

                for (std::size_t ip = 0; ip < phis.size(); ++ip) {
                    const EvalPhi ev = evaluate_phi(
                        q0,
                        phis[ip],
                        s0,
                        inv_lambda,
                        ewald_tol,
                        request.detector
                    );
                    if (!ev.ok) {
                        continue;
                    }
                    if (n_on == 0) {
                        phi_first = phis[ip];
                        x_mid = ev.x_mm;
                        y_mid = ev.y_mm;
                    }
                    phi_last = phis[ip];
                    // Refresh mid coords toward centre of observed interval.
                    if (n_on > 0) {
                        x_mid = 0.5 * (x_mid + ev.x_mm);
                        y_mid = 0.5 * (y_mid + ev.y_mm);
                    }
                    if (ip == 0) {
                        first_on = true;
                    }
                    if (ip + 1 == phis.size()) {
                        last_on = true;
                    }
                    ++n_on;
                }
                if (n_on == 0) {
                    continue;
                }

                PredictedReflection spot;
                spot.h = h;
                spot.k = k;
                spot.l = l;
                spot.x_mm = x_mid;
                spot.y_mm = y_mid;
                spot.phi_begin_deg = phi_first;
                spot.phi_end_deg = phi_last;
                if (still) {
                    spot.full = true;
                } else {
                    // Partial if rocking touches either oscillation endpoint.
                    spot.full = !(first_on || last_on);
                }
                spots.push_back(spot);
            }
        }
    }

    std::sort(spots.begin(), spots.end(), [](const auto& a, const auto& b) {
        if (a.h != b.h) {
            return a.h < b.h;
        }
        if (a.k != b.k) {
            return a.k < b.k;
        }
        return a.l < b.l;
    });
    PredictionResult result;
    result.reflections = std::move(spots);
    return result;
}

void ReflectionPredictor::write_dump(
    const std::filesystem::path& path,
    const std::vector<PredictedReflection>& spots
) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error(
            "ReflectionPredictor::write_dump: cannot open " + path.string()
        );
    }
    out << "PREDICT_DUMP n=" << spots.size() << '\n';
    out.setf(std::ios::fixed);
    for (const auto& s : spots) {
        out.precision(6);
        out << s.h << ' ' << s.k << ' ' << s.l << ' ' << s.x_mm << ' ' << s.y_mm
            << ' ';
        out.precision(4);
        out << s.phi_begin_deg << ' ' << s.phi_end_deg << ' '
            << (s.full ? 1 : 0) << '\n';
    }
}

std::vector<PredictedReflection> ReflectionPredictor::read_dump(
    const std::filesystem::path& path
) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error(
            "ReflectionPredictor::read_dump: cannot open " + path.string()
        );
    }
    std::string header;
    if (!std::getline(in, header)) {
        throw std::runtime_error("ReflectionPredictor::read_dump: empty file");
    }
    if (header.rfind("PREDICT_DUMP", 0) != 0) {
        throw std::runtime_error(
            "ReflectionPredictor::read_dump: missing PREDICT_DUMP header"
        );
    }

    std::vector<PredictedReflection> spots;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        std::istringstream iss(line);
        PredictedReflection s;
        int full_flag = 1;
        if (!(iss >> s.h >> s.k >> s.l >> s.x_mm >> s.y_mm >> s.phi_begin_deg >>
              s.phi_end_deg >> full_flag)) {
            throw std::runtime_error(
                "ReflectionPredictor::read_dump: bad line: " + line
            );
        }
        s.full = (full_flag != 0);
        spots.push_back(s);
    }
    return spots;
}

}  // namespace mosflm::diffraction
