#!/usr/bin/env python3
"""Independent stills Ewald-sphere predictor used to freeze the PREDICT oracle.

This is intentionally **not** the C++ implementation: it is a pure-Python
reference so ``testdata/oracle/predict_still.dump`` can be regenerated and the
C++ ``ReflectionPredictor`` differentially tested against it.

Full ``ipmosflm`` PREDICT capture still requires a known-indexed image set
(see ``testdata/oracle/README.md``). Until then this analytical dump is the
regression oracle for the stills vertical slice.

Usage::

    python3 tools/predict_still_oracle.py > testdata/oracle/predict_still.dump
"""

from __future__ import annotations

import math
import sys
from typing import List, Tuple


# --- Frozen fixture (thaumatin-like + simple detector) ---
A, B, C = 58.0, 58.0, 150.0
ALPHA = BETA = GAMMA = 90.0
WAVELENGTH = 1.0  # Å
DISTANCE_MM = 250.0
BEAM_X_MM = 0.0
BEAM_Y_MM = 0.0
PIXEL_X = PIXEL_Y = 0.172  # mm
NFAST, NSLOW = 1475, 1679
PHI_DEG = 0.0
DMIN = 3.5  # Å — keep dump small for CTest
EWALD_TOL = 0.0025  # Å⁻¹


def deg2rad(d: float) -> float:
    return d * math.pi / 180.0


def orth_matrix(
    a: float, b: float, c: float, alpha: float, beta: float, gamma: float
) -> List[List[float]]:
    """Fractional → Cartesian (Å); a along x, b in xy. Matches mosflm::UnitCell."""
    d2r = 1.0 / 57.29578
    ca, cb, cg = math.cos(alpha * d2r), math.cos(beta * d2r), math.cos(gamma * d2r)
    sg = math.sin(gamma * d2r)
    arg = 1.0 + 2.0 * ca * cb * cg - ca * ca - cb * cb - cg * cg
    V = a * b * c * math.sqrt(arg)
    return [
        [a, b * cg, c * cb],
        [0.0, b * sg, c * (ca - cb * cg) / sg],
        [0.0, 0.0, V / (a * b * sg)],
    ]


def mat_inv_3(M: List[List[float]]) -> List[List[float]]:
    a, b, c = M[0]
    d, e, f = M[1]
    g, h, i = M[2]
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    if abs(det) < 1e-18:
        raise RuntimeError("singular")
    inv = [
        [(e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det],
        [(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det],
        [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det],
    ]
    return inv


def mat_T(M: List[List[float]]) -> List[List[float]]:
    return [[M[j][i] for j in range(3)] for i in range(3)]


def mat_vec(M: List[List[float]], v: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (
        M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2],
        M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2],
        M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2],
    )


def rot_x(phi_deg: float, v: Tuple[float, float, float]) -> Tuple[float, float, float]:
    c = math.cos(deg2rad(phi_deg))
    s = math.sin(deg2rad(phi_deg))
    x, y, z = v
    return (x, c * y - s * z, s * y + c * z)


def d_spacing(h: int, k: int, l: int, Astar: List[List[float]]) -> float:
    q = mat_vec(Astar, (float(h), float(k), float(l)))
    q2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2]
    if q2 <= 0.0:
        return 1e30
    return 1.0 / math.sqrt(q2)


def lab_to_detector(
    s1: Tuple[float, float, float], distance: float, bx: float, by: float
) -> Tuple[bool, float, float]:
    if abs(s1[2]) < 1e-15:
        return False, 0.0, 0.0
    t = distance / s1[2]
    if t <= 0.0:
        return False, 0.0, 0.0
    return True, s1[0] * t + bx, s1[1] * t + by


def on_detector(x_mm: float, y_mm: float) -> bool:
    xf = x_mm / PIXEL_X
    yf = y_mm / PIXEL_Y
    # Origin at beam centre (0,0); accept if within half-detector about beam.
    return abs(xf) < 0.5 * NFAST and abs(yf) < 0.5 * NSLOW


def main() -> int:
    A_mat = orth_matrix(A, B, C, ALPHA, BETA, GAMMA)
    Astar = mat_T(mat_inv_3(A_mat))  # UB with R = I
    s0 = (0.0, 0.0, 1.0 / WAVELENGTH)
    inv_lambda = 1.0 / WAVELENGTH

    max_h = int(math.ceil(A / DMIN)) + 1
    max_k = int(math.ceil(B / DMIN)) + 1
    max_l = int(math.ceil(C / DMIN)) + 1

    spots = []
    for h in range(-max_h, max_h + 1):
        for k in range(-max_k, max_k + 1):
            for l in range(-max_l, max_l + 1):
                if h == 0 and k == 0 and l == 0:
                    continue
                if d_spacing(h, k, l, Astar) < DMIN:
                    continue
                q0 = mat_vec(Astar, (float(h), float(k), float(l)))
                q = rot_x(PHI_DEG, q0)
                s1 = (s0[0] + q[0], s0[1] + q[1], s0[2] + q[2])
                s1_len = math.sqrt(s1[0] ** 2 + s1[1] ** 2 + s1[2] ** 2)
                if abs(s1_len - inv_lambda) > EWALD_TOL:
                    continue
                ok, x_mm, y_mm = lab_to_detector(s1, DISTANCE_MM, BEAM_X_MM, BEAM_Y_MM)
                if not ok or not on_detector(x_mm, y_mm):
                    continue
                spots.append((h, k, l, x_mm, y_mm, PHI_DEG, PHI_DEG, 1))

    spots.sort(key=lambda t: (t[0], t[1], t[2]))
    print("PREDICT_DUMP n=%d" % len(spots))
    for h, k, l, x, y, pb, pe, full in spots:
        print("%d %d %d %.6f %.6f %.4f %.4f %d" % (h, k, l, x, y, pb, pe, full))
    return 0


if __name__ == "__main__":
    sys.exit(main())
