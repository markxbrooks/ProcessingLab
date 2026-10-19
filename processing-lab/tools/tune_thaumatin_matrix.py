#!/usr/bin/env python3
"""Short-window MOSFLM MATRIX fine-tune vs DIALS mean-I CC.

Holds CELL/BEAM/DISTANCE/PIXEL fixed. Searches small lab-frame missets
about the starting AMAT and a few axis-sign variants. Writes the winner
back into thaumatin_process.in when --write is set.

Usage (from ProcessingLab root)::

  python3 processing-lab/tools/tune_thaumatin_matrix.py --write
"""
from __future__ import annotations

import argparse
import math
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np

try:
    import gemmi
except ImportError as exc:  # pragma: no cover
    raise SystemExit("gemmi required (ccp4-python or pip install gemmi)") from exc

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_IN = ROOT / "thaumatin_process.in"
DEFAULT_PL = ROOT / "bin" / "ProcessingLab"
DEFAULT_IMAGE = Path.home() / "Downloads" / "thaumatin_9_1_master.h5"
DEFAULT_DIALS = (
    Path.home()
    / "Downloads"
    / "thaumatin_9_1_master"
    / "2026-10-15-20-29-46-dials"
    / "integrated.mtz"
)

HEADER = """\
! Auto-generated short-window MATRIX tune
SYMMETRY 75
WAVELENGTH 0.97949913
DISTANCE 199.78169135
PIXEL 0.075 0.075
BEAM 166.15771234 172.97150533
MOSAIC 0.15
SEPARATION 0.75 0.75
PEAKHALF 5
DMIN 1.15
DMAX 200.0
MINFRAC 0.5
CELL 57.93557545 57.93557545 150.27023041 90.0 90.0 90.0
MATRIX
{m00:.14g} {m01:.14g} {m02:.14g}
{m10:.14g} {m11:.14g} {m12:.14g}
{m20:.14g} {m21:.14g} {m22:.14g}
PROCESS {first} TO {last}
GO
EXIT
"""


def parse_amat(text: str) -> np.ndarray:
    m = re.search(
        r"MATRIX\s*\n"
        r"([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\n"
        r"([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\n"
        r"([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)",
        text,
    )
    if not m:
        raise ValueError("MATRIX block not found")
    vals = [float(m.group(i)) for i in range(1, 10)]
    return np.array(vals, dtype=float).reshape(3, 3)


def rot_axis(axis: str, deg: float) -> np.ndarray:
    a = math.radians(deg)
    c, s = math.cos(a), math.sin(a)
    if axis == "x":
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=float)
    if axis == "y":
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype=float)
    if axis == "z":
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=float)
    raise ValueError(axis)


def sign_variants(A: np.ndarray) -> list[tuple[str, np.ndarray]]:
    """Axis-sign / permutation-ish probes used historically."""
    out = [("base", A.copy())]
    Sz = np.diag([1.0, 1.0, -1.0])
    Sy = np.diag([1.0, -1.0, 1.0])
    Sx = np.diag([-1.0, 1.0, 1.0])
    out.append(("Sz_A", Sz @ A))
    out.append(("Sy_A", Sy @ A))
    out.append(("Sx_A", Sx @ A))
    out.append(("A_Sz", A @ Sz))
    return out


def mean_I_by_hkl(mtz: gemmi.Mtz, batch_lo: int, batch_hi: int) -> dict:
    h = mtz.column_with_label("H")
    k = mtz.column_with_label("K")
    l = mtz.column_with_label("L")
    inten = mtz.column_with_label("I")
    batch = mtz.column_with_label("BATCH")
    acc: dict[tuple[int, int, int], list[float]] = {}
    for i in range(mtz.nreflections):
        b = int(round(batch[i]))
        if b < batch_lo or b > batch_hi:
            continue
        key = (int(round(h[i])), int(round(k[i])), int(round(l[i])))
        acc.setdefault(key, []).append(float(inten[i]))
    return {k: float(np.mean(v)) for k, v in acc.items()}


def pearson_cc(a: dict, b: dict) -> tuple[float, int]:
    keys = sorted(set(a) & set(b))
    if len(keys) < 20:
        return float("nan"), len(keys)
    xa = np.array([a[k] for k in keys], dtype=float)
    xb = np.array([b[k] for k in keys], dtype=float)
    if np.std(xa) < 1e-12 or np.std(xb) < 1e-12:
        return float("nan"), len(keys)
    return float(np.corrcoef(xa, xb)[0, 1]), len(keys)


def median_rot_by_batch(mtz: gemmi.Mtz, batches: list[int]) -> dict[int, float]:
    batch = mtz.column_with_label("BATCH")
    rot = mtz.column_with_label("ROT")
    out = {}
    bb = np.array([batch[i] for i in range(mtz.nreflections)])
    rr = np.array([rot[i] for i in range(mtz.nreflections)])
    for bn in batches:
        m = bb == bn
        if m.any():
            out[bn] = float(np.median(rr[m]))
    return out


def run_process(
    pl: Path,
    image: Path,
    amat: np.ndarray,
    first: int,
    last: int,
    out_mtz: Path,
) -> None:
    script = HEADER.format(
        m00=amat[0, 0],
        m01=amat[0, 1],
        m02=amat[0, 2],
        m10=amat[1, 0],
        m11=amat[1, 1],
        m12=amat[1, 2],
        m20=amat[2, 0],
        m21=amat[2, 1],
        m22=amat[2, 2],
        first=first,
        last=last,
    )
    with tempfile.NamedTemporaryFile("w", suffix=".in", delete=False) as fh:
        fh.write(script)
        in_path = Path(fh.name)
    try:
        cmd = [
            str(pl),
            "HKLOUT",
            str(out_mtz),
            "--image",
            str(image),
            str(in_path),
        ]
        subprocess.run(
            cmd,
            check=True,
            cwd=str(ROOT),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    finally:
        in_path.unlink(missing_ok=True)


def write_amat_to_in(path: Path, amat: np.ndarray) -> None:
    text = path.read_text()
    block = (
        "MATRIX\n"
        f"{amat[0, 0]:.14g} {amat[0, 1]:.14g} {amat[0, 2]:.14g}\n"
        f"{amat[1, 0]:.14g} {amat[1, 1]:.14g} {amat[1, 2]:.14g}\n"
        f"{amat[2, 0]:.14g} {amat[2, 1]:.14g} {amat[2, 2]:.14g}\n"
    )
    new, n = re.subn(
        r"MATRIX\s*\n"
        r"[-\d.eE+]+\s+[-\d.eE+]+\s+[-\d.eE+]+\s*\n"
        r"[-\d.eE+]+\s+[-\d.eE+]+\s+[-\d.eE+]+\s*\n"
        r"[-\d.eE+]+\s+[-\d.eE+]+\s+[-\d.eE+]+\s*\n",
        block,
        text,
        count=1,
    )
    if n != 1:
        raise RuntimeError(f"failed to rewrite MATRIX in {path}")
    path.write_text(new)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_path", type=Path, default=DEFAULT_IN)
    ap.add_argument("--pl", type=Path, default=DEFAULT_PL)
    ap.add_argument("--image", type=Path, default=DEFAULT_IMAGE)
    ap.add_argument("--dials", type=Path, default=DEFAULT_DIALS)
    ap.add_argument("--first", type=int, default=2)
    ap.add_argument("--last", type=int, default=26)
    ap.add_argument("--step", type=float, default=0.1)
    ap.add_argument("--span", type=float, default=0.5)
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--check-rot", action="store_true")
    args = ap.parse_args()

    base = parse_amat(args.in_path.read_text())
    dials = gemmi.read_mtz_file(str(args.dials))
    dials_mean = mean_I_by_hkl(dials, args.first, args.last)

    if args.check_rot:
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "rot.mtz"
            run_process(args.pl, args.image, base, args.first, args.last, out)
            pl = gemmi.read_mtz_file(str(out))
            batches = list(range(args.first, min(args.first + 3, args.last + 1)))
            pr = median_rot_by_batch(pl, batches)
            dr = median_rot_by_batch(dials, batches)
            print("ROT median check (PL vs DIALS):")
            for bn in batches:
                if bn in pr and bn in dr:
                    d = abs(pr[bn] - dr[bn])
                    print(f"  batch {bn}: PL={pr[bn]:.4f} DI={dr[bn]:.4f} |d|={d:.4f}")
                    if d > 0.05:
                        print("  WARNING: |dROT| > 0.05°")
        return 0

    def score(name: str, A: np.ndarray, out: Path) -> tuple[float, int]:
        run_process(args.pl, args.image, A, args.first, args.last, out)
        pl = gemmi.read_mtz_file(str(out))
        pl_mean = mean_I_by_hkl(pl, args.first, args.last)
        return pearson_cc(pl_mean, dials_mean)

    # Prefer high overlap: reject candidates with n < min_frac * max_n seen.
    min_frac = 0.5
    results: list[tuple[str, np.ndarray, float, int]] = []
    work = Path(tempfile.mkdtemp(prefix="tune_amat_"))
    print(f"Phase 1: sign variants on frames {args.first}-{args.last}")
    print(f"work dir: {work}")

    for i, (name, A) in enumerate(sign_variants(base)):
        out = work / f"s{i:02d}.mtz"
        try:
            cc, n = score(name, A, out)
        except Exception as exc:  # noqa: BLE001
            print(f"  FAIL {name}: {exc}")
            continue
        tag = f"cc={cc:.4f} n={n}" if cc == cc else f"cc=nan n={n}"
        print(f"  {name}: {tag}")
        if cc == cc:
            results.append((name, A, cc, n))

    max_n = max((n for _, _, _, n in results), default=0)
    n_floor = int(min_frac * max_n)
    print(f"n_floor={n_floor} (min_frac={min_frac} of max_n={max_n})")

    def pick_best(
        pool: list[tuple[str, np.ndarray, float, int]]
    ) -> tuple[str, np.ndarray, float, int]:
        ok = [r for r in pool if r[3] >= n_floor]
        if not ok:
            ok = pool
        return max(ok, key=lambda r: r[2])

    best_name, best_A, best_cc, best_n = pick_best(results)
    print(f"Phase-1 pick: {best_name} cc={best_cc:.4f} n={best_n}")

    angles: list[float] = []
    a = -args.span
    while a <= args.span + 1e-12:
        if abs(a) >= 1e-12:
            angles.append(round(a, 6))
        a += args.step

    print(f"Phase 2: misset about {best_name} (±{args.span}° step {args.step}°)")
    seed_A = best_A.copy()
    seed_name = best_name
    idx = 0
    for ax in ("x", "y", "z"):
        for ang in angles:
            name = f"{seed_name}+R{ax}({ang})"
            A = rot_axis(ax, ang) @ seed_A
            out = work / f"m{idx:03d}.mtz"
            idx += 1
            try:
                cc, n = score(name, A, out)
            except Exception as exc:  # noqa: BLE001
                print(f"  FAIL {name}: {exc}")
                continue
            if cc != cc:
                continue
            results.append((name, A, cc, n))
            if n >= n_floor and cc > best_cc:
                best_cc, best_name, best_A, best_n = cc, name, A, n
                print(f"  NEW BEST {name}: cc={cc:.4f} n={n}")

    # Combined Rx,Ry,Rz at best 1D angles (small product search).
    print("Phase 3: joint ±step about current best")
    seed_A = best_A.copy()
    seed_name = best_name
    for dx in (-args.step, 0.0, args.step):
        for dy in (-args.step, 0.0, args.step):
            for dz in (-args.step, 0.0, args.step):
                if abs(dx) + abs(dy) + abs(dz) < 1e-12:
                    continue
                name = f"{seed_name}+Rxyz({dx},{dy},{dz})"
                A = rot_axis("z", dz) @ rot_axis("y", dy) @ rot_axis("x", dx) @ seed_A
                out = work / f"j{idx:03d}.mtz"
                idx += 1
                try:
                    cc, n = score(name, A, out)
                except Exception as exc:  # noqa: BLE001
                    print(f"  FAIL {name}: {exc}")
                    continue
                if cc != cc:
                    continue
                results.append((name, A, cc, n))
                if n >= n_floor and cc > best_cc:
                    best_cc, best_name, best_A, best_n = cc, name, A, n
                    print(f"  NEW BEST {name}: cc={cc:.4f} n={n}")

    print(f"BEST {best_name}: mean-I CC={best_cc:.4f} on {best_n} common HKL")
    print("AMAT:")
    for r in range(3):
        print(f"  {best_A[r,0]:.14g} {best_A[r,1]:.14g} {best_A[r,2]:.14g}")

    if args.write and best_cc == best_cc:
        write_amat_to_in(args.in_path, best_A)
        print(f"Wrote MATRIX to {args.in_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
