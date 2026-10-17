/**
 * C ABI for Phase C call-outs from Fortran / C.
 *
 * Stable POD layout so ``iso_c_binding`` or a thin ``*_`` Fortran wrapper can
 * invoke the C++ stills predictor and summation integrator without pulling
 * C++ types across the boundary.
 *
 * Build: linked into ``mosflm_geometry`` (C++ TU). Callers need only this
 * header and the library.
 */
#ifndef MOSFLM_C_API_H
#define MOSFLM_C_API_H

#ifdef __cplusplus
extern "C" {
#endif

/** Maximum reflections returned by a single ``mosflm_predict_still`` call. */
#define MOSFLM_C_MAX_SPOTS 8192

typedef struct mosflm_cell6 {
    double a;
    double b;
    double c;
    double alpha_deg;
    double beta_deg;
    double gamma_deg;
} mosflm_cell6;

typedef struct mosflm_detector_c {
    double distance_mm;
    double beam_x_mm;
    double beam_y_mm;
    double pixel_size_x_mm;
    double pixel_size_y_mm;
    int nfast;
    int nslow;
    double wavelength_A;
} mosflm_detector_c;

typedef struct mosflm_predicted_c {
    int h;
    int k;
    int l;
    double x_mm;
    double y_mm;
    double phi_begin_deg;
    double phi_end_deg;
    int full; /* 1 = full, 0 = partial */
} mosflm_predicted_c;

typedef struct mosflm_integrated_c {
    int h;
    int k;
    int l;
    double x_mm;
    double y_mm;
    double intensity;
    double sigma;
    int ok;
} mosflm_integrated_c;

/**
 * Stills PREDICT with identity orientation (UB from cell, R = I).
 *
 * \param cell        Direct cell
 * \param det         Detector geometry (beam-relative mm frame)
 * \param phi_deg     Still φ
 * \param dmin_A      Resolution limit
 * \param ewald_tol   Ewald tolerance (Å⁻¹); use 0.0025 for analytical fixture
 * \param out         Caller-allocated array (≥ MOSFLM_C_MAX_SPOTS recommended)
 * \param out_cap     Capacity of ``out``
 * \param out_n       Number of spots written
 * \return 0 on success, non-zero on error
 */
int mosflm_predict_still(
    const mosflm_cell6* cell,
    const mosflm_detector_c* det,
    double phi_deg,
    double dmin_A,
    double ewald_tol,
    mosflm_predicted_c* out,
    int out_cap,
    int* out_n
);

/**
 * Apply analytical stills fixture defaults into ``cell`` / ``det``
 * (thaumatin-like cell, 1475×1679 detector) for smoke tests.
 */
void mosflm_analytical_still_fixture(
    mosflm_cell6* cell,
    mosflm_detector_c* det,
    double* phi_deg,
    double* dmin_A,
    double* ewald_tol
);

/**
 * Lorentz–polarization factor for a spot at detector mm (beam-centred).
 * Unpolarized beam: ``P = (1 + cos² 2θ) / 2``, ``L = 1 / sin(2θ)`` (clamped).
 * Returns ``L * P`` (dimensionless scale for intensities).
 */
double mosflm_lp_factor(
    const mosflm_detector_c* det,
    double x_mm,
    double y_mm
);

/**
 * Box summation integrate over a raster (row-major, slow×fast = ``nslow×nfast``).
 *
 * \param apply_lp  Non-zero → scale intensity/sigma by ``mosflm_lp_factor``
 * \return 0 on success, non-zero on error (1 = bad args, 2 = truncated)
 */
int mosflm_integrate_summation(
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
);

#ifdef __cplusplus
}
#endif

#endif /* MOSFLM_C_API_H */
