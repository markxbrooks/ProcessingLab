#include "mosflm/legacy/MosflmCApi.h"

#include "MosflmCApiInternal.hpp"

#include <cmath>
#include <vector>

namespace {

void pack_and_predict_still(
    float cell[6],
    double* distance_mm,
    double* wavelength_A,
    double* pixel_x_mm,
    double* pixel_y_mm,
    int* nfast,
    int* nslow,
    double* beam_x_mm,
    double* beam_y_mm,
    double* phi_deg,
    double* dmin_A,
    double* ewald_tol,
    const mosflm_ub9* ub,
    int* h,
    int* k,
    int* l,
    double* x_mm,
    double* y_mm,
    double* phi_beg,
    double* phi_end,
    int* ifull,
    const int* nmax,
    int* nspot,
    int* ierr
) {
    if (ierr) {
        *ierr = MOSFLM_ERROR_INVALID_ARGUMENT;
    }
    if (nspot) {
        *nspot = 0;
    }
    if (!cell || !distance_mm || !wavelength_A || !pixel_x_mm || !pixel_y_mm ||
        !nfast || !nslow || !beam_x_mm || !beam_y_mm || !phi_deg || !dmin_A ||
        !ewald_tol || !h || !k || !l || !x_mm || !y_mm || !phi_beg || !phi_end ||
        !ifull || !nmax || !nspot || !ierr || *nmax <= 0) {
        return;
    }

    try {
        *ierr = MOSFLM_SUCCESS;
        *nspot = 0;

        mosflm_cell6 c{};
        c.a = cell[0];
        c.b = cell[1];
        c.c = cell[2];
        c.alpha_deg = cell[3];
        c.beta_deg = cell[4];
        c.gamma_deg = cell[5];

        mosflm_detector_c d{};
        d.distance_mm = *distance_mm;
        d.wavelength_A = *wavelength_A;
        d.pixel_size_x_mm = *pixel_x_mm;
        d.pixel_size_y_mm = *pixel_y_mm;
        d.nfast = *nfast;
        d.nslow = *nslow;
        d.beam_x_mm = *beam_x_mm;
        d.beam_y_mm = *beam_y_mm;

        std::vector<mosflm_predicted_c> buf(static_cast<std::size_t>(*nmax));
        const int rc = mosflm_predict_still(
            &c, &d, ub, nullptr, *phi_deg, *dmin_A, *ewald_tol, buf.data(),
            *nmax, nspot
        );
        *ierr = rc;
        for (int i = 0; i < *nspot; ++i) {
            h[i] = buf[static_cast<std::size_t>(i)].h;
            k[i] = buf[static_cast<std::size_t>(i)].k;
            l[i] = buf[static_cast<std::size_t>(i)].l;
            x_mm[i] = buf[static_cast<std::size_t>(i)].x_mm;
            y_mm[i] = buf[static_cast<std::size_t>(i)].y_mm;
            phi_beg[i] = buf[static_cast<std::size_t>(i)].phi_begin_deg;
            phi_end[i] = buf[static_cast<std::size_t>(i)].phi_end_deg;
            ifull[i] = buf[static_cast<std::size_t>(i)].full;
        }
    } catch (...) {
        *nspot = 0;
        *ierr = MOSFLM_ERROR_INTERNAL;
    }
}

}  // namespace

extern "C" void predict_still_cpp_(
    float cell[6],
    double* distance_mm,
    double* wavelength_A,
    double* pixel_x_mm,
    double* pixel_y_mm,
    int* nfast,
    int* nslow,
    double* beam_x_mm,
    double* beam_y_mm,
    double* phi_deg,
    double* dmin_A,
    double* ewald_tol,
    int* h,
    int* k,
    int* l,
    double* x_mm,
    double* y_mm,
    double* phi_beg,
    double* phi_end,
    int* ifull,
    const int* nmax,
    int* nspot,
    int* ierr
) {
    pack_and_predict_still(
        cell, distance_mm, wavelength_A, pixel_x_mm, pixel_y_mm, nfast, nslow,
        beam_x_mm, beam_y_mm, phi_deg, dmin_A, ewald_tol, nullptr, h, k, l,
        x_mm, y_mm, phi_beg, phi_end, ifull, nmax, nspot, ierr
    );
}

extern "C" void predict_still_psi_cpp_(
    float cell[6],
    double* distance_mm,
    double* wavelength_A,
    double* pixel_x_mm,
    double* pixel_y_mm,
    int* nfast,
    int* nslow,
    double* beam_x_mm,
    double* beam_y_mm,
    double* phi_deg,
    double* dmin_A,
    double* ewald_tol,
    double psi_deg[3],
    int* h,
    int* k,
    int* l,
    double* x_mm,
    double* y_mm,
    double* phi_beg,
    double* phi_end,
    int* ifull,
    const int* nmax,
    int* nspot,
    int* ierr
) {
    if (ierr) {
        *ierr = MOSFLM_ERROR_INVALID_ARGUMENT;
    }
    if (nspot) {
        *nspot = 0;
    }
    if (!cell || !psi_deg) {
        return;
    }

    mosflm_cell6 c{};
    c.a = cell[0];
    c.b = cell[1];
    c.c = cell[2];
    c.alpha_deg = cell[3];
    c.beta_deg = cell[4];
    c.gamma_deg = cell[5];

    mosflm_ub9 ub{};
    const int ub_rc = mosflm_ub_from_psi(&c, psi_deg, &ub);
    if (ub_rc != MOSFLM_SUCCESS) {
        if (ierr) {
            *ierr = ub_rc;
        }
        return;
    }

    pack_and_predict_still(
        cell, distance_mm, wavelength_A, pixel_x_mm, pixel_y_mm, nfast, nslow,
        beam_x_mm, beam_y_mm, phi_deg, dmin_A, ewald_tol, &ub, h, k, l, x_mm,
        y_mm, phi_beg, phi_end, ifull, nmax, nspot, ierr
    );
}

extern "C" void predict_osc_amat_cpp_(
    float cell[6],
    double* distance_mm,
    double* wavelength_A,
    double* pixel_x_mm,
    double* pixel_y_mm,
    int* nfast,
    int* nslow,
    double* beam_x_mm,
    double* beam_y_mm,
    double* phi_start_deg,
    double* phi_end_deg,
    double* dmin_A,
    double* ewald_tol,
    double* phi_step_deg,
    float amat_colmaj[9],
    int* h,
    int* k,
    int* l,
    double* x_mm,
    double* y_mm,
    double* phi_beg,
    double* phi_end,
    int* ifull,
    const int* nmax,
    int* nspot,
    int* ierr
) {
    if (ierr) {
        *ierr = MOSFLM_ERROR_INVALID_ARGUMENT;
    }
    if (nspot) {
        *nspot = 0;
    }
    if (!cell || !distance_mm || !wavelength_A || !pixel_x_mm || !pixel_y_mm ||
        !nfast || !nslow || !beam_x_mm || !beam_y_mm || !phi_start_deg ||
        !phi_end_deg || !dmin_A || !ewald_tol || !amat_colmaj || !h || !k ||
        !l || !x_mm || !y_mm || !phi_beg || !phi_end || !ifull || !nmax ||
        !nspot || !ierr || *nmax <= 0) {
        return;
    }

    mosflm_ub9 ub{};
    const int ub_rc = mosflm_ub_from_mosflm_amat(amat_colmaj, &ub);
    if (ub_rc != MOSFLM_SUCCESS) {
        *ierr = ub_rc;
        return;
    }

    const double step =
        (phi_step_deg && std::isfinite(*phi_step_deg) && *phi_step_deg > 0.0)
            ? *phi_step_deg
            : 0.05;

    try {
        *ierr = MOSFLM_SUCCESS;
        *nspot = 0;

        mosflm_cell6 c{};
        c.a = cell[0];
        c.b = cell[1];
        c.c = cell[2];
        c.alpha_deg = cell[3];
        c.beta_deg = cell[4];
        c.gamma_deg = cell[5];

        mosflm_detector_c d{};
        d.distance_mm = *distance_mm;
        d.wavelength_A = *wavelength_A;
        d.pixel_size_x_mm = *pixel_x_mm;
        d.pixel_size_y_mm = *pixel_y_mm;
        d.nfast = *nfast;
        d.nslow = *nslow;
        d.beam_x_mm = *beam_x_mm;
        d.beam_y_mm = *beam_y_mm;

        std::vector<mosflm_predicted_c> buf(static_cast<std::size_t>(*nmax));
        const int rc = mosflm_predict_osc(
            &c, &d, &ub, nullptr, *phi_start_deg, *phi_end_deg, *dmin_A,
            *ewald_tol, step, buf.data(), *nmax, nspot
        );
        *ierr = rc;
        for (int i = 0; i < *nspot; ++i) {
            h[i] = buf[static_cast<std::size_t>(i)].h;
            k[i] = buf[static_cast<std::size_t>(i)].k;
            l[i] = buf[static_cast<std::size_t>(i)].l;
            x_mm[i] = buf[static_cast<std::size_t>(i)].x_mm;
            y_mm[i] = buf[static_cast<std::size_t>(i)].y_mm;
            phi_beg[i] = buf[static_cast<std::size_t>(i)].phi_begin_deg;
            phi_end[i] = buf[static_cast<std::size_t>(i)].phi_end_deg;
            ifull[i] = buf[static_cast<std::size_t>(i)].full;
        }
    } catch (...) {
        *nspot = 0;
        *ierr = MOSFLM_ERROR_INTERNAL;
    }
}
