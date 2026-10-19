#pragma once

#include "mosflm/legacy/MosflmCApi.h"

#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/geometry/UnitCell.hpp"

namespace mosflm::legacy::capi {

inline void zero_out_n(int* out_n) {
    if (out_n) {
        *out_n = 0;
    }
}

[[nodiscard]] bool finite6(const mosflm_cell6& c);

[[nodiscard]] bool finite_detector(const mosflm_detector_c& d);

[[nodiscard]] bool finite_scalar(double x);

[[nodiscard]] bool valid_predict_buffers(mosflm_predicted_c* out, int out_cap);

[[nodiscard]] bool valid_integrate_dims(int nfast, int nslow);

[[nodiscard]] bool dims_match_detector(
    const mosflm_detector_c& det,
    int nfast,
    int nslow
);

[[nodiscard]] mosflm::UnitCell to_cell(const mosflm_cell6& c);

[[nodiscard]] mosflm::DetectorGeometry to_det(const mosflm_detector_c& d);

[[nodiscard]] mosflm::OrientationMatrix to_orientation(
    const mosflm::UnitCell& cell,
    const mosflm_ub9* ub
);

[[nodiscard]] const char* spacegroup_or_default(const char* spacegroup);

/**
 * Core predict. Caller has already zeroed ``out_n`` and validated pointers
 * sufficiently for count-only vs write modes.
 */
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
);

void fill_cell_det_from_analytical_still(
    mosflm_cell6* cell,
    mosflm_detector_c* det,
    double* phi_deg,
    double* dmin_A,
    double* ewald_tol
);

void fill_cell_det_from_analytical_osc(
    mosflm_cell6* cell,
    mosflm_detector_c* det,
    double* phi_start_deg,
    double* phi_end_deg,
    double* dmin_A,
    double* ewald_tol,
    double* phi_step_deg
);

}  // namespace mosflm::legacy::capi
