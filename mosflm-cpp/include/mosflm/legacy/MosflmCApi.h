/**
 * C ABI for Phase C call-outs from Fortran / C.
 *
 * Stable POD layout so ``iso_c_binding`` or a thin F77 wrapper can invoke the
 * C++ stills/oscillation predictor and summation integrator without pulling
 * C++ types across the boundary.
 *
 * Contracts
 * ---------
 * - Status codes: see ``MOSFLM_SUCCESS`` / ``MOSFLM_ERROR_*`` below.
 * - ``out_n`` is always the **number of elements written** into ``out``.
 *   When ``out_n`` is non-null it is set to 0 before validation on every
 *   exported predict/integrate entry (including failure paths).
 * - Truncation: if more spots/reflections are produced than ``out_cap``,
 *   the first ``out_cap`` are written, ``out_n == out_cap``, and the return
 *   is ``MOSFLM_ERROR_BUFFER_TOO_SMALL``.
 * - Count-only: ``out_cap == 0`` and ``out == NULL`` is allowed for predict;
 *   ``out_n`` receives the total count and the return is ``MOSFLM_SUCCESS``.
 * - Exceptions never escape ``extern "C"`` entry points (mapped to
 *   ``MOSFLM_ERROR_INTERNAL``).
 * - Pixel layout for integrate: row-major, index ``row * nfast + col``
 *   (fast contiguous). Beam ``x/y`` mm are beam-relative.
 * - Orientation: pass ``mosflm_ub9`` (lab UB, ``r* = UB * h``, row-major).
 *   ``NULL`` selects the canonical fixture orientation
 *   ``from_cell_and_rotation(cell, I)``. Spindle φ is applied about +X on
 *   top of UB; the unit cell alone does not define crystal orientation.
 * - Space group: ``NULL`` or empty ⇒ ``"P 1"``. The predictor does **not**
 *   apply systematic absences yet; SG is metadata / future filtering.
 * - F77 symbol ``predict_still_cpp_`` assumes **gfortran** trailing-underscore
 *   calling convention.
 *
 * Build: linked into ``mosflm_geometry`` (C++ TU). Callers need only this
 * header and the library.
 */
#ifndef MOSFLM_C_API_H
#define MOSFLM_C_API_H

#ifdef __cplusplus
extern "C" {
#endif

/** Maximum reflections returned by a single predict call (recommended cap). */
#define MOSFLM_C_MAX_SPOTS 8192

/** Maximum allowed image dimension (fast or slow) for integrate. */
#define MOSFLM_C_MAX_IMAGE_DIM 10000

enum {
    MOSFLM_SUCCESS = 0,
    MOSFLM_ERROR_INVALID_ARGUMENT = 1,
    MOSFLM_ERROR_BUFFER_TOO_SMALL = 2,
    MOSFLM_ERROR_INTERNAL = 3
};

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

/** Lab UB matrix, row-major 3×3 (``r* = UB * h``). */
typedef struct mosflm_ub9 {
    double m[9];
} mosflm_ub9;

/**
 * Build lab UB from cell + MOSFLM missetting angles PSI (degrees).
 *
 * Uses ``UB = R(PSI) · A*`` with ``R = Rz · Ry · Rx`` (``ROTMAT``).
 * ``psi_deg`` is length 3: phix, phiy, phiz.
 *
 * \return ``MOSFLM_SUCCESS`` or ``MOSFLM_ERROR_INVALID_ARGUMENT``
 */
int mosflm_ub_from_psi(
    const mosflm_cell6* cell,
    const double psi_deg[3],
    mosflm_ub9* out
);

/**
 * Copy a row-major 3×3 matrix (MOSFLM ``MATRIX`` / A-matrix style) into
 * ``out``. ``m`` must be length 9.
 *
 * \return ``MOSFLM_SUCCESS`` or ``MOSFLM_ERROR_INVALID_ARGUMENT``
 */
int mosflm_ub_from_matrix9(const double m[9], mosflm_ub9* out);

/**
 * Build C++ predictor UB from a MOSFLM Orientation Matrix [A]
 * (Fortran column-major ``AMAT(3,3)`` as 9 contiguous floats) via
 * :func:`mosflm_amat_to_cpp_ub`.
 *
 * \return ``MOSFLM_SUCCESS`` or ``MOSFLM_ERROR_INVALID_ARGUMENT``
 */
int mosflm_ub_from_mosflm_amat(const float amat_colmaj[9], mosflm_ub9* out);

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
 * Stills PREDICT.
 *
 * \param ub          Lab UB or NULL for canonical cell orientation (R = I)
 * \param spacegroup  Hermann–Mauguin name or NULL/empty for ``"P 1"``
 * \param out_n       Number written (or total in count-only mode)
 * \return ``MOSFLM_*`` status
 */
int mosflm_predict_still(
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
 * Analytical oscillation fixture (``predict_osc.dump``): φ ∈ [0, 0.5]°,
 * dmin = 5 Å. Writes ``phi_start`` / ``phi_end``.
 */
void mosflm_analytical_osc_fixture(
    mosflm_cell6* cell,
    mosflm_detector_c* det,
    double* phi_start_deg,
    double* phi_end_deg,
    double* dmin_A,
    double* ewald_tol,
    double* phi_step_deg
);

/**
 * Oscillation PREDICT. ``phi_end == phi_start`` is a still.
 */
int mosflm_predict_osc(
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
);

/**
 * F77-callable stills PREDICT (gfortran trailing underscore).
 *
 * Production path only: caller supplies ``CELL(6)`` and detector scalars.
 * Use ``mosflm_analytical_still_fixture`` from C/tests to fill inputs first.
 */
void predict_still_cpp_(
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
);

/**
 * F77 stills PREDICT with missetting PSI (degrees): builds UB via
 * ``mosflm_ub_from_psi`` then calls the same path as ``mosflm_predict_still``.
 * ``psi_deg`` is length 3 (phix, phiy, phiz). Zero PSI matches
 * ``predict_still_cpp_`` (canonical R = I).
 */
void predict_still_psi_cpp_(
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
);

/**
 * Oscillation PREDICT with MOSFLM AMAT (column-major 3×3) mapped to C++ UB.
 * ``phi_end == phi_start`` is a still. ``phi_step_deg`` may be NULL → 0.05°.
 */
void predict_osc_amat_cpp_(
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
);

/**
 * Lorentz–polarization factor (delegates to domain ``LpCorrection``).
 * Unpolarized: ``P = (1 + cos² 2θ) / 2``, ``L = 1 / |sin 2θ|``; near-beam
 * clamp returns 1.0 (not physically valid on-axis).
 */
double mosflm_lp_factor(
    const mosflm_detector_c* det,
    double x_mm,
    double y_mm
);

/**
 * Box summation integrate (row-major ``nslow×nfast``).
 *
 * \param apply_lp  Non-zero → scale intensity/sigma inside the integrator
 * \param n_preds   May be 0 with ``preds == NULL``
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
