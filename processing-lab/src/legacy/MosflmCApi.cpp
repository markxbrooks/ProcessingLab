#include "mosflm/legacy/MosflmCApi.h"

#include "MosflmCApiInternal.hpp"

#include "mosflm/diffraction/LpCorrection.hpp"
#include "mosflm/diffraction/SummationIntegrator.hpp"
#include "mosflm/geometry/LabFrame.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/geometry/Rotmat.hpp"
#include "mosflm/geometry/UnitCell.hpp"
#include "mosflm/image/DiffractionImage.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace {

bool finite_triple(const double v[3]) {
    return v && std::isfinite(v[0]) && std::isfinite(v[1]) &&
           std::isfinite(v[2]);
}

void fill_ub9(const Eigen::Matrix3d& M, mosflm_ub9* out) {
    out->m[0] = M(0, 0);
    out->m[1] = M(0, 1);
    out->m[2] = M(0, 2);
    out->m[3] = M(1, 0);
    out->m[4] = M(1, 1);
    out->m[5] = M(1, 2);
    out->m[6] = M(2, 0);
    out->m[7] = M(2, 1);
    out->m[8] = M(2, 2);
}

}  // namespace

extern "C" int mosflm_ub_from_psi(
    const mosflm_cell6* cell,
    const double psi_deg[3],
    mosflm_ub9* out
) {
    if (!cell || !out || !finite_triple(psi_deg) ||
        !mosflm::legacy::capi::finite6(*cell)) {
        return MOSFLM_ERROR_INVALID_ARGUMENT;
    }
    try {
        const mosflm::UnitCell uc{
            cell->a, cell->b, cell->c, cell->alpha_deg, cell->beta_deg,
            cell->gamma_deg};
        const Eigen::Matrix3d R = mosflm::rotmat_from_missetting(
            psi_deg[0], psi_deg[1], psi_deg[2], true
        );
        const auto om =
            mosflm::OrientationMatrix::from_cell_and_rotation(uc, R);
        fill_ub9(om.ub(), out);
        return MOSFLM_SUCCESS;
    } catch (...) {
        return MOSFLM_ERROR_INTERNAL;
    }
}

extern "C" int mosflm_ub_from_matrix9(const double m[9], mosflm_ub9* out) {
    if (!m || !out) {
        return MOSFLM_ERROR_INVALID_ARGUMENT;
    }
    for (int i = 0; i < 9; ++i) {
        if (!std::isfinite(m[i])) {
            return MOSFLM_ERROR_INVALID_ARGUMENT;
        }
        out->m[i] = m[i];
    }
    return MOSFLM_SUCCESS;
}

extern "C" int mosflm_ub_from_mosflm_amat(
    const float amat_colmaj[9],
    mosflm_ub9* out
) {
    if (!amat_colmaj || !out) {
        return MOSFLM_ERROR_INVALID_ARGUMENT;
    }
    for (int i = 0; i < 9; ++i) {
        if (!std::isfinite(amat_colmaj[i])) {
            return MOSFLM_ERROR_INVALID_ARGUMENT;
        }
    }
    try {
        Eigen::Matrix3d A;
        A << amat_colmaj[0], amat_colmaj[3], amat_colmaj[6],
            amat_colmaj[1], amat_colmaj[4], amat_colmaj[7],
            amat_colmaj[2], amat_colmaj[5], amat_colmaj[8];
        fill_ub9(mosflm::mosflm_amat_to_cpp_ub(A), out);
        return MOSFLM_SUCCESS;
    } catch (...) {
        return MOSFLM_ERROR_INTERNAL;
    }
}

extern "C" void mosflm_analytical_still_fixture(
    mosflm_cell6* cell,
    mosflm_detector_c* det,
    double* phi_deg,
    double* dmin_A,
    double* ewald_tol
) {
    try {
        mosflm::legacy::capi::fill_cell_det_from_analytical_still(
            cell, det, phi_deg, dmin_A, ewald_tol
        );
    } catch (...) {
        // Fixtures are best-effort; leave outputs untouched on failure.
    }
}

extern "C" void mosflm_analytical_osc_fixture(
    mosflm_cell6* cell,
    mosflm_detector_c* det,
    double* phi_start_deg,
    double* phi_end_deg,
    double* dmin_A,
    double* ewald_tol,
    double* phi_step_deg
) {
    try {
        mosflm::legacy::capi::fill_cell_det_from_analytical_osc(
            cell, det, phi_start_deg, phi_end_deg, dmin_A, ewald_tol,
            phi_step_deg
        );
    } catch (...) {
    }
}

extern "C" int mosflm_predict_still(
    const mosflm_cell6* cell,
    const mosflm_detector_c* det,
    const mosflm_ub9* ub,
    const char* spacegroup,
    double phi_deg,
    double dmin_A,
    double ewald_tol,
    mosflm_predicted_c* out,
    int out_cap,
    int* out_n
) {
    mosflm::legacy::capi::zero_out_n(out_n);
    if (!out_n || !cell || !det ||
        !mosflm::legacy::capi::valid_predict_buffers(out, out_cap)) {
        return MOSFLM_ERROR_INVALID_ARGUMENT;
    }
    try {
        return mosflm::legacy::capi::run_predict(
            *cell, *det, ub, spacegroup, phi_deg, phi_deg, dmin_A, ewald_tol,
            0.0, out, out_cap, out_n
        );
    } catch (...) {
        *out_n = 0;
        return MOSFLM_ERROR_INTERNAL;
    }
}

extern "C" int mosflm_predict_osc(
    const mosflm_cell6* cell,
    const mosflm_detector_c* det,
    const mosflm_ub9* ub,
    const char* spacegroup,
    double phi_start_deg,
    double phi_end_deg,
    double dmin_A,
    double ewald_tol,
    double phi_step_deg,
    mosflm_predicted_c* out,
    int out_cap,
    int* out_n
) {
    mosflm::legacy::capi::zero_out_n(out_n);
    if (!out_n || !cell || !det ||
        !mosflm::legacy::capi::valid_predict_buffers(out, out_cap)) {
        return MOSFLM_ERROR_INVALID_ARGUMENT;
    }
    try {
        return mosflm::legacy::capi::run_predict(
            *cell, *det, ub, spacegroup, phi_start_deg, phi_end_deg, dmin_A,
            ewald_tol, phi_step_deg, out, out_cap, out_n
        );
    } catch (...) {
        *out_n = 0;
        return MOSFLM_ERROR_INTERNAL;
    }
}

extern "C" double mosflm_lp_factor(
    const mosflm_detector_c* det,
    double x_mm,
    double y_mm
) {
    try {
        if (!det || !mosflm::legacy::capi::finite_detector(*det)) {
            return 1.0;
        }
        if (!mosflm::legacy::capi::finite_scalar(x_mm) ||
            !mosflm::legacy::capi::finite_scalar(y_mm)) {
            return 1.0;
        }
        return mosflm::diffraction::lp_factor(
            mosflm::legacy::capi::to_det(*det), x_mm, y_mm
        );
    } catch (...) {
        return 1.0;
    }
}

extern "C" int mosflm_integrate_summation(
    const mosflm_detector_c* det,
    const int* pixels,
    int nfast,
    int nslow,
    const mosflm_predicted_c* preds,
    int n_preds,
    int box_half_px,
    int bg_frame_px,
    int apply_lp,
    mosflm_integrated_c* out,
    int out_cap,
    int* out_n
) {
    mosflm::legacy::capi::zero_out_n(out_n);
    if (!out_n || !det || !pixels || !out || out_cap <= 0 || n_preds < 0 ||
        (n_preds > 0 && !preds) || box_half_px < 0 || bg_frame_px < 0 ||
        !mosflm::legacy::capi::valid_integrate_dims(nfast, nslow) ||
        !mosflm::legacy::capi::finite_detector(*det) ||
        !mosflm::legacy::capi::dims_match_detector(*det, nfast, nslow)) {
        return MOSFLM_ERROR_INVALID_ARGUMENT;
    }
    try {
        mosflm::image::DiffractionImage image;
        image.metadata.width = nfast;
        image.metadata.height = nslow;
        const std::size_t npix =
            static_cast<std::size_t>(nfast) * static_cast<std::size_t>(nslow);
        image.pixels.resize(npix);
        for (std::size_t i = 0; i < npix; ++i) {
            image.pixels[i] = static_cast<std::int32_t>(pixels[i]);
        }

        mosflm::DetectorGeometry g = mosflm::legacy::capi::to_det(*det);
        g.nfast = nfast;
        g.nslow = nslow;

        std::vector<mosflm::diffraction::PredictedReflection> pred_v;
        pred_v.reserve(static_cast<std::size_t>(n_preds));
        for (int i = 0; i < n_preds; ++i) {
            mosflm::diffraction::PredictedReflection p;
            p.h = preds[i].h;
            p.k = preds[i].k;
            p.l = preds[i].l;
            p.x_mm = preds[i].x_mm;
            p.y_mm = preds[i].y_mm;
            p.phi_begin_deg = preds[i].phi_begin_deg;
            p.phi_end_deg = preds[i].phi_end_deg;
            p.full = preds[i].full != 0;
            pred_v.push_back(p);
        }

        mosflm::diffraction::IntegrationSettings settings;
        settings.box_half_px = box_half_px;
        settings.bg_frame_px = bg_frame_px;
        settings.apply_lp = apply_lp != 0;
        const mosflm::diffraction::IntegrationRequest ireq{
            image, g, pred_v, settings};

        const auto integrated =
            mosflm::diffraction::SummationIntegrator{}.integrate(ireq);
        const int n = static_cast<int>(
            std::min(integrated.size(), static_cast<std::size_t>(out_cap))
        );
        for (int i = 0; i < n; ++i) {
            const auto& r = integrated[static_cast<std::size_t>(i)];
            out[i].h = r.h;
            out[i].k = r.k;
            out[i].l = r.l;
            out[i].x_mm = r.x_mm;
            out[i].y_mm = r.y_mm;
            out[i].intensity = r.intensity;
            out[i].sigma = r.sigma;
            out[i].ok = r.ok ? 1 : 0;
        }
        *out_n = n;
        return (integrated.size() > static_cast<std::size_t>(out_cap))
                   ? MOSFLM_ERROR_BUFFER_TOO_SMALL
                   : MOSFLM_SUCCESS;
    } catch (...) {
        *out_n = 0;
        return MOSFLM_ERROR_INTERNAL;
    }
}
