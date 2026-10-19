#!/usr/bin/env python3
"""Independent oscillation Ewald predictor → ``predict_osc.dump``.

Matches ``ReflectionPredictor`` oscillation rules (φ about +X, sample step,
full = not touching oscillation endpoints).

Usage::

    python3 tools/predict_osc_oracle.py > testdata/oracle/predict_osc.dump
"""

from __future__ import annotations

import math
import sys
from typing import List, Tuple

A, B, C = 58.0, 58.0, 150.0
ALPHA = BETA = GAMMA = 90.0
WAVELENGTH = 1.0
DISTANCE_MM = 250.0
BEAM_X_MM = BEAM_Y_MM = 0.0
PIXEL_X = PIXEL_Y = 0.172
NFAST, NSLOW = 1475, 1679
PHI_START, PHI_END = 0.0, 0.5
PHI_STEP = 0.05
DMIN = 5.0
EWALD_TOL = 0.0025


def deg2rad(d: float) -> float:
    return d * math.pi / 180.0


def orth_matrix(a, b, c, alpha, beta, gamma):
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


def mat_inv_3(M):
    a, b, c = M[0]
    d, e, f = M[1]
    g, h, i = M[2]
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    return [
        [(e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det],
        [(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det],
        [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det],
    ]


def mat_T(M):
    return [[M[j][i] for j in range(3)] for i in range(3)]


def mat_vec(M, v):
    return (
        M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2],
        M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2],
        M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2],
    )


def rot_x(phi_deg, v):
    c = math.cos(deg2rad(phi_deg))
    s = math.sin(deg2rad(phi_deg))
    x, y, z = v
    return (x, c * y - s * z, s * y + c * z)


def d_spacing(h, k, l, Astar):
    q = mat_vec(Astar, (float(h), float(k), float(l)))
    q2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2]
    return 1e30 if q2 <= 0.0 else 1.0 / math.sqrt(q2)


def lab_to_detector(s1, distance, bx, by):
    if abs(s1[2]) < 1e-15:
        return False, 0.0, 0.0
    t = distance / s1[2]
    if t <= 0.0:
        return False, 0.0, 0.0
    return True, s1[0] * t + bx, s1[1] * t + by


def on_detector(x_mm, y_mm):
    return abs(x_mm / PIXEL_X) < 0.5 * NFAST and abs(y_mm / PIXEL_Y) < 0.5 * NSLOW


def main() -> int:
    A_mat = orth_matrix(A, B, C, ALPHA, BETA, GAMMA)
    Astar = mat_T(mat_inv_3(A_mat))
    s0 = (0.0, 0.0, 1.0 / WAVELENGTH)
    inv_lambda = 1.0 / WAVELENGTH

    phis = []
    phi = PHI_START
    while phi <= PHI_END + 1e-9:
        phis.append(min(phi, PHI_END))
        phi += PHI_STEP
    if abs(phis[-1] - PHI_END) > 1e-9:
        phis.append(PHI_END)

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
                n_on = 0
                phi_first = phi_last = 0.0
                x_mid = y_mid = 0.0
                first_on = last_on = False
                for ip, phi in enumerate(phis):
                    q = rot_x(phi, q0)
                    s1 = (s0[0] + q[0], s0[1] + q[1], s0[2] + q[2])
                    s1_len = math.sqrt(s1[0] ** 2 + s1[1] ** 2 + s1[2] ** 2)
                    if abs(s1_len - inv_lambda) > EWALD_TOL:
                        continue
                    ok, x_mm, y_mm = lab_to_detector(
                        s1, DISTANCE_MM, BEAM_X_MM, BEAM_Y_MM
                    )
                    if not ok or not on_detector(x_mm, y_mm):
                        continue
                    if n_on == 0:
                        phi_first = phi
                        x_mid, y_mid = x_mm, y_mm
                    else:
                        x_mid = 0.5 * (x_mid + x_mm)
                        y_mid = 0.5 * (y_mid + y_mm)
                    phi_last = phi
                    if ip == 0:
                        first_on = True
                    if ip + 1 == len(phis):
                        last_on = True
                    n_on += 1
                if n_on == 0:
                    continue
                full = 0 if (first_on or last_on) else 1
                spots.append(
                    (h, k, l, x_mid, y_mid, phi_first, phi_last, full)
                )

    spots.sort(key=lambda t: (t[0], t[1], t[2]))
    print("PREDICT_DUMP n=%d" % len(spots))
    for h, k, l, x, y, pb, pe, full in spots:
        print("%d %d %d %.6f %.6f %.4f %.4f %d" % (h, k, l, x, y, pb, pe, full))
    return 0


if __name__ == "__main__":
    sys.exit(main())
