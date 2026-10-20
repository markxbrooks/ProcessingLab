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
[[nodiscard]] Eigen::Matrix3d mosflm_amat_to_cpp_ub(
    const Eigen::Matrix3d& A_mosflm,
    double wavelength_A = 1.0
);

/// Inverse of :func:`mosflm_amat_to_cpp_ub` (CONTROL ``MATRIX`` snapshot).
[[nodiscard]] Eigen::Matrix3d cpp_ub_to_mosflm_amat(
    const Eigen::Matrix3d& ub_cpp,
    double wavelength_A = 1.0
);

/**
 * DPSCHOOSE / ``spots_to_rlp_xyz`` frame → C++ predictor frame.
 *
 * ``spots_to_rlp_xyz`` builds RLPs in the MOSFLM REFIX frame (beam +z,
 * spindle +y, detector rotated by the ``.spt`` OMEGA, fast axis mirrored when
 * INVERTX). The predictor rotates about +x with the detector mapped through
 * ``DetectorGeometry::lab_to_detector_mm`` (CCOM). Exact identity, for every
 * ω and φ, between the inverse of the RLP build and the predictor:
 *
 * ``r_cpp = Q · r_xyz`` with ``CCOM = ccom_from_spt_omega(ω, invertx)``,
 * ``Q = [[0,1,0],[-1,0,0],[0,0,1]]`` for INVERTX (MOSFLM/ADSC ``.spt``) and
 * ``Q = [[0,-1,0],[-1,0,0],[0,0,1]]`` for non-inverted files (ProcessingLab
 * FINDSPOTS). Verified numerically on thaumatin (ω = 270, INVERTX 0).
 */
[[nodiscard]] Eigen::Matrix3d dps_xyz_to_cpp_lab(bool invertx);

/**
 * MOSFLM ``CCOM`` (``DetectorGeometry::ccom_deg``) for a ``.spt`` / NeXus
 * detector omega so PROCESS predicts in the frame AUTOINDEX indexed in:
 * ``−ω`` when INVERTX, ``ω + 180`` otherwise; normalised to [0, 360).
 */
[[nodiscard]] double ccom_from_spt_omega(double omega_deg, bool invertx = true);

/**
 * Convert a reciprocal basis in the ``spots_to_rlp_xyz`` frame
 * (``om_t_xyz`` columns = a*, b*, c* in Å⁻¹, i.e. ``h = om_t_xyz⁻¹ · r_xyz``)
 * to a MOSFLM-convention A-matrix for CONTROL / PROCESS / XPARM.
 */
[[nodiscard]] Eigen::Matrix3d xyz_reciprocal_to_mosflm_amat(
    const Eigen::Matrix3d& om_t_xyz,
    double wavelength_A,
    bool invertx
);

/// Inverse of :func:`xyz_reciprocal_to_mosflm_amat`.
[[nodiscard]] Eigen::Matrix3d mosflm_amat_to_xyz_reciprocal(
    const Eigen::Matrix3d& A_mosflm,
    double wavelength_A,
    bool invertx
);

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

/**
 * Same as :func:`cpp_mm_to_mosflm_xdet`, applying MOSFLM ``YSCALE`` to the
 * slow-axis pixel pitch (``pixel_y_eff = pixel_y / yscale``). Matches
 * MeasurementBox / StandardProfile integration centres.
 */
void cpp_mm_to_mosflm_xdet_yscale(
    double x_mm,
    double y_mm,
    double beam_x_mm,
    double beam_y_mm,
    double pixel_size_x_mm,
    double pixel_size_y_mm,
    double yscale,
    double& xdet_px,
    double& ydet_px
);

}  // namespace mosflm
