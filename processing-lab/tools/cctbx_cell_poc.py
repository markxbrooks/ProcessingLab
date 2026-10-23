#!/usr/bin/env python
"""CCTBX selective-eval POC: unit cell + change-of-basis vs MOSFLM / Gemmi.

Prefer Phenix ``cctbx.python`` (Python 2.7 in Phenix 1.20)::

    /Applications/phenix-*/build/bin/cctbx.python tools/cctbx_cell_poc.py

See ``docs/CCTBX_EVAL.md`` for adopt/defer decisions. This script is a
standalone evaluator; it is **not** linked into ``processinglab_geometry``.
"""

from __future__ import print_function

import math
import sys


def mosflm_volume(a, b, c, alpha, beta, gamma):
    """Match ``mosflm::UnitCell::volume`` (MOSFLM deg2rad = 1/57.29578)."""
    d2r = 1.0 / 57.29578
    ca = math.cos(alpha * d2r)
    cb = math.cos(beta * d2r)
    cc = math.cos(gamma * d2r)
    arg = 1.0 + 2.0 * ca * cb * cc - ca * ca - cb * cb - cc * cc
    if arg <= 0.0:
        return 0.0
    return a * b * c * math.sqrt(arg)


def mosflm_reciprocal(a, b, c, alpha, beta, gamma):
    """Match ``mosflm::UnitCell::reciprocal`` / INVCEL."""
    d2r = 1.0 / 57.29578
    r2d = 57.29578
    ca = [
        math.cos(alpha * d2r),
        math.cos(beta * d2r),
        math.cos(gamma * d2r),
    ]
    sa = [
        math.sin(alpha * d2r),
        math.sin(beta * d2r),
        math.sin(gamma * d2r),
    ]
    cell = [a, b, c]
    V = mosflm_volume(a, b, c, alpha, beta, gamma)
    out = [0.0] * 6
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        ri = cell[j] * cell[k] * sa[i] / V
        arg = r2d * math.acos((ca[j] * ca[k] - ca[i]) / (sa[j] * sa[k]))
        rounded = int(round(arg))
        if abs(arg - rounded) < 0.0009:
            arg = float(rounded)
        out[i] = ri
        out[i + 3] = arg
    return tuple(out)


def near(a, b, tol):
    return abs(a - b) <= tol


def main():
    try:
        from cctbx import uctbx, sgtbx
        from scitbx import matrix
    except ImportError as exc:
        print(
            "FAIL: cctbx not importable (%s). Run with Phenix cctbx.python "
            "or install cctbx." % exc,
            file=sys.stderr,
        )
        return 2

    # Orthorhombic cell used across processing-lab geometry tests.
    params = (58.0, 72.0, 90.0, 90.0, 90.0, 90.0)
    uc = uctbx.unit_cell(params)
    m_vol = mosflm_volume(*params)
    m_rec = mosflm_reciprocal(*params)
    c_rec = uc.reciprocal_parameters()

    print("=== Cell comparison (direct) ===")
    print("params     :", params)
    print("cctbx vol  :", uc.volume())
    print("mosflm vol :", m_vol)
    print("cctbx recip:", c_rec)
    print("mosflm recip:", m_rec)

    ok = True
    if not near(uc.volume(), m_vol, 1e-4):
        print("FAIL: volume mismatch")
        ok = False
    for i, (ca, ma) in enumerate(zip(c_rec, m_rec)):
        tol = 1e-8 if i < 3 else 1e-4
        if not near(ca, ma, tol):
            print("FAIL: reciprocal param[%d] cctbx=%s mosflm=%s" % (i, ca, ma))
            ok = False

    # Gemmi agreement is covered by CTest gemmi_crystal; restate expected a*.
    expected_astar = 1.0 / 58.0
    if not near(c_rec[0], expected_astar, 1e-10):
        print("FAIL: a* vs 1/a")
        ok = False
    else:
        print("Gemmi/MOSFLM a* check (1/a): OK (%.10f)" % expected_astar)

    # DPS lattice-character transform #2 (hR):
    # {1,-1,0; -1,0,1; -1,-1,-1} from bravais_transform.c tm[2].
    print("=== Change-of-basis (DPS character 2 matrix) ===")
    m = matrix.sqr((1, -1, 0, -1, 0, 1, -1, -1, -1))
    print("matrix det:", m.determinant())
    try:
        cb = sgtbx.change_of_basis_op(sgtbx.rt_mx(m, (0, 0, 0)))
    except Exception:
        cb = sgtbx.change_of_basis_op("a-b, -a+c, -a-b-c")
    print("cb as_abc :", cb.as_abc())

    cubic = uctbx.unit_cell((50.0, 50.0, 50.0, 90.0, 90.0, 90.0))
    new_uc = cubic.change_basis(cb)
    print("cubic     :", cubic.parameters())
    print("transformed:", new_uc.parameters())
    print(
        "new volume :",
        new_uc.volume(),
        "vs |det|*V =",
        abs(m.determinant()) * cubic.volume(),
    )

    if not near(new_uc.volume(), abs(m.determinant()) * cubic.volume(), 1e-6):
        print("FAIL: transformed volume vs |det| * V")
        ok = False

    id_cb = sgtbx.change_of_basis_op("a,b,c")
    same = uc.change_basis(id_cb)
    if same.parameters() != params:
        print("FAIL: identity CoB changed cell", same.parameters())
        ok = False
    else:
        print("identity CoB: OK")

    if ok:
        print("CCTBX_POC_OK")
        return 0
    print("CCTBX_POC_FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
