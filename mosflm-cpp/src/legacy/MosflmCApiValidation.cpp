#include "MosflmCApiInternal.hpp"

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace mosflm::legacy::capi {

bool finite6(const mosflm_cell6& c) {
    return std::isfinite(c.a) && std::isfinite(c.b) && std::isfinite(c.c) &&
           std::isfinite(c.alpha_deg) && std::isfinite(c.beta_deg) &&
           std::isfinite(c.gamma_deg) && c.a > 0.0 && c.b > 0.0 && c.c > 0.0;
}

bool finite_detector(const mosflm_detector_c& d) {
    return std::isfinite(d.distance_mm) && d.distance_mm > 0.0 &&
           std::isfinite(d.beam_x_mm) && std::isfinite(d.beam_y_mm) &&
           std::isfinite(d.pixel_size_x_mm) && d.pixel_size_x_mm > 0.0 &&
           std::isfinite(d.pixel_size_y_mm) && d.pixel_size_y_mm > 0.0 &&
           std::isfinite(d.wavelength_A) && d.wavelength_A > 0.0;
}

bool finite_scalar(double x) {
    return std::isfinite(x);
}

bool valid_predict_buffers(mosflm_predicted_c* out, int out_cap) {
    if (out_cap < 0) {
        return false;
    }
    if (out_cap == 0) {
        return out == nullptr;  // count-only
    }
    return out != nullptr;
}

bool valid_integrate_dims(int nfast, int nslow) {
    if (nfast <= 0 || nslow <= 0) {
        return false;
    }
    if (nfast > MOSFLM_C_MAX_IMAGE_DIM || nslow > MOSFLM_C_MAX_IMAGE_DIM) {
        return false;
    }
    const std::uint64_t nf = static_cast<std::uint64_t>(nfast);
    const std::uint64_t ns = static_cast<std::uint64_t>(nslow);
    if (nf * ns > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
        return false;
    }
    return true;
}

bool dims_match_detector(
    const mosflm_detector_c& det,
    int nfast,
    int nslow
) {
    if (det.nfast != 0 && det.nfast != nfast) {
        return false;
    }
    if (det.nslow != 0 && det.nslow != nslow) {
        return false;
    }
    return true;
}

mosflm::UnitCell to_cell(const mosflm_cell6& c) {
    return mosflm::UnitCell{
        c.a, c.b, c.c, c.alpha_deg, c.beta_deg, c.gamma_deg
    };
}

mosflm::DetectorGeometry to_det(const mosflm_detector_c& d) {
    mosflm::DetectorGeometry g;
    g.distance_mm = d.distance_mm;
    g.beam_x_mm = d.beam_x_mm;
    g.beam_y_mm = d.beam_y_mm;
    g.pixel_size_x_mm = d.pixel_size_x_mm;
    g.pixel_size_y_mm = d.pixel_size_y_mm;
    g.nfast = d.nfast;
    g.nslow = d.nslow;
    g.wavelength_A = d.wavelength_A;
    return g;
}

mosflm::OrientationMatrix to_orientation(
    const mosflm::UnitCell& cell,
    const mosflm_ub9* ub
) {
    if (!ub) {
        return mosflm::OrientationMatrix::from_cell_and_rotation(cell);
    }
    Eigen::Matrix3d M;
    M << ub->m[0], ub->m[1], ub->m[2], ub->m[3], ub->m[4], ub->m[5], ub->m[6],
        ub->m[7], ub->m[8];
    return mosflm::OrientationMatrix(M);
}

const char* spacegroup_or_default(const char* spacegroup) {
    if (!spacegroup || spacegroup[0] == '\0') {
        return "P 1";
    }
    return spacegroup;
}

int run_predict(
    const mosflm_cell6& cell,
    const mosflm_detector_c& det,
    const mosflm_ub9* ub,
    const char* spacegroup,
    double phi_start,
    double phi_end,
    double dmin_A,
    double ewald_tol,
    double phi_step,
    mosflm_predicted_c* out,
    int out_cap,
    int* out_n
) {
    *out_n = 0;
    if (!finite6(cell) || !finite_detector(det) || !finite_scalar(phi_start) ||
        !finite_scalar(phi_end) || !finite_scalar(dmin_A) ||
        !finite_scalar(ewald_tol) || !(dmin_A > 0.0) || phi_end < phi_start) {
        return MOSFLM_ERROR_INVALID_ARGUMENT;
    }
    if (ub) {
        for (int i = 0; i < 9; ++i) {
            if (!finite_scalar(ub->m[i])) {
                return MOSFLM_ERROR_INVALID_ARGUMENT;
            }
        }
    }
    if (phi_step != 0.0 && !finite_scalar(phi_step)) {
        return MOSFLM_ERROR_INVALID_ARGUMENT;
    }

    mosflm::diffraction::PredictionRequest req;
    const mosflm::UnitCell uc = to_cell(cell);
    req.crystal = mosflm::crystal::CrystalMetadata::from_mosflm_cell(
        uc, spacegroup_or_default(spacegroup)
    );
    req.orientation = to_orientation(uc, ub);
    req.detector = to_det(det);
    req.scan.phi_start_deg = phi_start;
    req.scan.phi_end_deg = phi_end;
    req.settings.dmin_A = dmin_A;
    req.settings.ewald_tolerance_Ainv =
        (ewald_tol > 0.0 && std::isfinite(ewald_tol)) ? ewald_tol : 0.0025;
    if (phi_step > 0.0) {
        req.settings.phi_step_deg = phi_step;
    }

    const auto result =
        mosflm::diffraction::ReflectionPredictor{}.predict(req);
    const auto& spots = result.reflections;
    if (out_cap == 0) {
        *out_n = static_cast<int>(spots.size());
        return MOSFLM_SUCCESS;
    }
    const int n = static_cast<int>(
        std::min(spots.size(), static_cast<std::size_t>(out_cap))
    );
    for (int i = 0; i < n; ++i) {
        const auto& s = spots[static_cast<std::size_t>(i)];
        out[i].h = s.h;
        out[i].k = s.k;
        out[i].l = s.l;
        out[i].x_mm = s.x_mm;
        out[i].y_mm = s.y_mm;
        out[i].phi_begin_deg = s.phi_begin_deg;
        out[i].phi_end_deg = s.phi_end_deg;
        out[i].full = s.full ? 1 : 0;
    }
    *out_n = n;
    return (spots.size() > static_cast<std::size_t>(out_cap))
               ? MOSFLM_ERROR_BUFFER_TOO_SMALL
               : MOSFLM_SUCCESS;
}

void fill_cell_det_from_analytical_still(
    mosflm_cell6* cell,
    mosflm_detector_c* det,
    double* phi_deg,
    double* dmin_A,
    double* ewald_tol
) {
    const auto req =
        mosflm::diffraction::ReflectionPredictor::analytical_still_fixture();
    if (cell) {
        cell->a = req.crystal.cell.a;
        cell->b = req.crystal.cell.b;
        cell->c = req.crystal.cell.c;
        cell->alpha_deg = req.crystal.cell.alpha;
        cell->beta_deg = req.crystal.cell.beta;
        cell->gamma_deg = req.crystal.cell.gamma;
    }
    if (det) {
        det->distance_mm = req.detector.distance_mm;
        det->beam_x_mm = req.detector.beam_x_mm;
        det->beam_y_mm = req.detector.beam_y_mm;
        det->pixel_size_x_mm = req.detector.pixel_size_x_mm;
        det->pixel_size_y_mm = req.detector.pixel_size_y_mm;
        det->nfast = req.detector.nfast;
        det->nslow = req.detector.nslow;
        det->wavelength_A = req.detector.wavelength_A;
    }
    if (phi_deg) {
        *phi_deg = req.scan.phi_start_deg;
    }
    if (dmin_A) {
        *dmin_A = req.settings.dmin_A;
    }
    if (ewald_tol) {
        *ewald_tol = req.settings.ewald_tolerance_Ainv;
    }
}

void fill_cell_det_from_analytical_osc(
    mosflm_cell6* cell,
    mosflm_detector_c* det,
    double* phi_start_deg,
    double* phi_end_deg,
    double* dmin_A,
    double* ewald_tol,
    double* phi_step_deg
) {
    const auto req = mosflm::diffraction::ReflectionPredictor::
        analytical_oscillation_fixture();
    if (cell) {
        cell->a = req.crystal.cell.a;
        cell->b = req.crystal.cell.b;
        cell->c = req.crystal.cell.c;
        cell->alpha_deg = req.crystal.cell.alpha;
        cell->beta_deg = req.crystal.cell.beta;
        cell->gamma_deg = req.crystal.cell.gamma;
    }
    if (det) {
        det->distance_mm = req.detector.distance_mm;
        det->beam_x_mm = req.detector.beam_x_mm;
        det->beam_y_mm = req.detector.beam_y_mm;
        det->pixel_size_x_mm = req.detector.pixel_size_x_mm;
        det->pixel_size_y_mm = req.detector.pixel_size_y_mm;
        det->nfast = req.detector.nfast;
        det->nslow = req.detector.nslow;
        det->wavelength_A = req.detector.wavelength_A;
    }
    if (phi_start_deg) {
        *phi_start_deg = req.scan.phi_start_deg;
    }
    if (phi_end_deg) {
        *phi_end_deg = req.scan.phi_end_deg;
    }
    if (dmin_A) {
        *dmin_A = req.settings.dmin_A;
    }
    if (ewald_tol) {
        *ewald_tol = req.settings.ewald_tolerance_Ainv;
    }
    if (phi_step_deg) {
        *phi_step_deg = req.settings.phi_step_deg;
    }
}

}  // namespace mosflm::legacy::capi
