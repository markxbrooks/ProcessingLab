#pragma once

#include <Eigen/Dense>

namespace mosflm {

/**
 * Laboratory-frame conventions.
 *
 * MOSFLM (``/REEKE/``): ``S0 ≈ (-1,0,0)`` is **antiparallel** to the X-ray beam
 * (photons travel along +X), spindle about +Z, detector downstream along +X.
 * Virtual detector ``XG`` ∥ lab Y, ``YG`` ∥ lab Z (see ``unoptimizable.f``
 * GENDATA comments).
 *
 * C++ :class:`ReflectionPredictor`: incident wavevector ``s0 = (0,0,1/λ)``
 * (photons along +Z), spindle about +X, detector at ``z = distance``.
 *
 * Right-handed map ``v_cpp = R * v_mosflm`` with
 * ``R = [[0,0,1],[0,-1,0],[1,0,0]]``:
 * MOSFLM +X → C++ +Z, +Y → −Y, +Z → +X.
 *
 * Empirically, applying ``R`` to the printed Orientation Matrix [A] recovers
 * every insulin MTZ HKL as a space-group / Friedel equivalent of a C++
 * prediction (exact Miller indices may differ by symmetry).
 */
[[nodiscard]] Eigen::Matrix3d mosflm_lab_to_cpp();

/**
 * Convert a MOSFLM orientation / A-matrix (rows = lab X,Y,Z; columns = A*,B*,C*)
 * into the C++ predictor UB (``r*_cpp = UB_cpp * h``).
 */
[[nodiscard]] Eigen::Matrix3d mosflm_amat_to_cpp_ub(const Eigen::Matrix3d& A_mosflm);

/**
 * Map C++ predictor detector millimetres (beam-relative, from
 * :meth:`DetectorGeometry::lab_to_detector_mm`) to MOSFLM ``XDET``/``YDET``
 * pixel coordinates (origin = image lower-left).
 *
 * GENDATA: ``XG`` ∥ lab Y, ``YG`` ∥ lab Z. After :func:`mosflm_lab_to_cpp`,
 * insulin known-indexed matching gives ``XG ↔ +y_mm``, ``YG ↔ +x_mm``
 * (XG sense is empirical; YG follows Z→X).
 *
 * ``beam_x_mm`` / ``beam_y_mm`` are the direct-beam position in the image
 * frame (SMV ``BEAM`` / MOSFLM ``XCEN``/``YCEN``), not the C++ predictor's
 * beam-centred origin.
 */
void cpp_mm_to_mosflm_xdet(
    double x_mm,
    double y_mm,
    double beam_x_mm,
    double beam_y_mm,
    double pixel_size_x_mm,
    double pixel_size_y_mm,
    double& xdet_px,
    double& ydet_px
);

}  // namespace mosflm
