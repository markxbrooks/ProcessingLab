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

}  // namespace mosflm
