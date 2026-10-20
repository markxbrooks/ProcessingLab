#!/usr/bin/env python3
"""Intensity benchmark: ProcessingLab PROCESS vs DIALS ``integrated.mtz``.

Runs the autoindex workflow (FINDSPOTS 1-50, AUTOINDEX DPS, PROCESS 1..N)
on the thaumatin Eiger master, then compares reflection by reflection with
the DIALS integration of the same images. Matching uses original (unreduced)
Miller indices and the nearest rotation angle within 2 degrees.

Usage (from ProcessingLab root)::

  python3 processing-lab/tools/bench_thaumatin_vs_dials.py            # 100 frames
  python3 processing-lab/tools/bench_thaumatin_vs_dials.py --frames 300
  python3 processing-lab/tools/bench_thaumatin_vs_dials.py --mtz process.mtz
"""
from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import gemmi
except ImportError as exc:  # pragma: no cover
    raise SystemExit("gemmi required (pip install gemmi)") from exc

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PL = ROOT / "bin" / "ProcessingLab"
DEFAULT_IMAGE = Path.home() / "Downloads" / "thaumatin_9_1_master.h5"
DEFAULT_DIALS = (
    Path.home()
    / "Downloads"
    / "thaumatin_9_1_master"
    / "2026-10-15-20-29-46-dials"
    / "integrated.mtz"
)
RES_BINS = [1.0, 1.6, 2.0, 2.5, 3.5, 5.0, 10.0, 1e9]
STRONG_ISIGI = 10.0


def run_processinglab(
    exe: Path, image: Path, frames: int, nproc: int, workdir: Path, extra: str
) -> Path:
    """Run the FINDSPOTS / AUTOINDEX / PROCESS workflow; return the MTZ path.

    :param exe: ProcessingLab executable.
    :param image: HDF5 master file.
    :param frames: Last frame to PROCESS (first is 1).
    :param nproc: OpenMP threads.
    :param workdir: Directory for the spot file, MTZ and log.
    :param extra: Extra CONTROL lines inserted before PROCESS.
    :return: Path of the written MTZ.
    """
    mtz = workdir / "bench.mtz"
    spt = workdir / "bench.spt"
    script = (
        f"NPROC {nproc}\n"
        f"FINDSPOTS THRESHOLD 5.0 MAXSPOTS 1500 FIND 1 FILE {spt} "
        "MULTIFRAME RANGE 1 50\nGO\n"
        "AUTOINDEX DPS THRESHOLD 15\nGO\n"
        f"{extra}"
        f"PROCESS 1 TO {frames}\nGO\nEXIT\n"
    )
    log = workdir / "bench.log"
    with log.open("w") as fh:
        subprocess.run(
            [str(exe), "--image", str(image), "HKLOUT", str(mtz)],
            input=script,
            text=True,
            stdout=fh,
            stderr=subprocess.STDOUT,
            check=True,
            cwd=workdir,
        )
    return mtz


def load(path: Path) -> pd.DataFrame:
    """Read an unmerged MTZ with original indices and resolution.

    :param path: MTZ file.
    :return: One row per observation.
    """
    mtz = gemmi.read_mtz_file(str(path))
    mtz.switch_to_original_hkl()
    df = pd.DataFrame(np.array(mtz, copy=False), columns=mtz.column_labels())
    df = df.astype({"H": int, "K": int, "L": int, "BATCH": int})
    df["d"] = mtz.cell.calculate_d_array(df[["H", "K", "L"]].to_numpy())
    return df


def cc(x: pd.Series, y: pd.Series) -> float:
    """Pearson correlation, NaN when fewer than 10 pairs.

    :param x: First series.
    :param y: Second series.
    :return: Correlation coefficient.
    """
    return float(np.corrcoef(x, y)[0, 1]) if len(x) >= 10 else float("nan")


def compare(pl: pd.DataFrame, di: pd.DataFrame) -> dict:
    """Compute the benchmark metrics.

    :param pl: ProcessingLab observations.
    :param di: DIALS observations.
    :return: Metric dictionary (JSON serialisable).
    """
    rot_hi = pl.ROT.max() + 0.2
    di = di[di.ROT <= rot_hi]
    # Only fully recorded ProcessingLab reflections (edges of the processed
    # range are partial by construction).
    n_partial = int((pl.FRACTIONCALC < 0.95).sum())
    pl = pl[pl.FRACTIONCALC >= 0.95]
    m = pd.merge_asof(
        pl.sort_values("ROT"),
        di.sort_values("ROT"),
        on="ROT",
        by=["H", "K", "L"],
        direction="nearest",
        tolerance=2.0,
        suffixes=("_pl", "_di"),
    ).dropna(subset=["I_di"])
    m["off_px"] = np.hypot(m.XDET_pl - m.XDET_di, m.YDET_pl - m.YDET_di)
    strong = m[m.I_di / m.SIGI_di > STRONG_ISIGI]
    full = pl.FRACTIONCALC >= 0.999

    # Reflection passages: same HKL, rows less than 3 degrees apart.
    s = pl.sort_values(["H", "K", "L", "ROT"])
    passage = (s.groupby(["H", "K", "L"]).ROT.diff().fillna(999.0) > 3.0).cumsum()
    multi_full = int(
        (s[full.loc[s.index]].groupby(passage[full.loc[s.index]]).size() > 1).sum()
    )

    def by_bins(frame: pd.DataFrame, key: pd.Series) -> dict:
        g = frame.assign(_k=key).groupby("_k", observed=True)
        return {
            str(k): {
                "n": int(len(v)),
                "ratio": float(np.median(v.I_pl / v.I_di)),
                "cc": cc(v.I_pl, v.I_di),
                "off_px": float(v.off_px.median()),
            }
            for k, v in g
        }

    return {
        "rows_pl": int(len(pl)),
        "rows_pl_partial_skipped": n_partial,
        "rows_dials_in_range": int(len(di)),
        "matched": int(len(m)),
        "strong_n": int(len(strong)),
        "cc_all": cc(m.I_pl, m.I_di),
        "cc_strong": cc(strong.I_pl, strong.I_di),
        "cc_strong_log": cc(
            np.log(np.abs(strong.I_pl) + 1.0), np.log(strong.I_di + 1.0)
        ),
        "median_ratio_strong": float(np.median(strong.I_pl / strong.I_di)),
        "median_offset_px": float(m.off_px.median()),
        "isigi_median_pl": float(np.median(pl.I / pl.SIGI.where(pl.SIGI > 0))),
        "isigi_median_dials": float(np.median(di.I / di.SIGI)),
        "zero_i_rows": int((pl.I == 0).sum()),
        "passages_multi_full": multi_full,
        "by_resolution": by_bins(strong, pd.cut(strong.d_pl, RES_BINS)),
        "by_rotation": by_bins(strong, (strong.ROT // 5 * 5).astype(int)),
    }


def print_report(r: dict) -> None:
    """Print the metrics as plain tables.

    :param r: Output of :func:`compare`.
    """
    for k in (
        "rows_pl", "rows_pl_partial_skipped", "rows_dials_in_range", "matched", "strong_n", "cc_all",
        "cc_strong", "cc_strong_log", "median_ratio_strong", "median_offset_px",
        "isigi_median_pl", "isigi_median_dials", "zero_i_rows",
        "passages_multi_full",
    ):
        v = r[k]
        print(f"  {k:22s} {v:.4f}" if isinstance(v, float) else f"  {k:22s} {v}")
    for title in ("by_resolution", "by_rotation"):
        print(f"\n  {title} (strong):  n  ratio  cc  offset_px")
        for key, v in r[title].items():
            print(
                f"    {key:>16s} {v['n']:7d} {v['ratio']:7.3f} "
                f"{v['cc']:6.3f} {v['off_px']:6.2f}"
            )


def main() -> None:
    """Command-line entry point."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--frames", type=int, default=100)
    ap.add_argument("--nproc", type=int, default=8)
    ap.add_argument("--exe", type=Path, default=DEFAULT_PL)
    ap.add_argument("--image", type=Path, default=DEFAULT_IMAGE)
    ap.add_argument("--dials", type=Path, default=DEFAULT_DIALS)
    ap.add_argument("--mtz", type=Path, help="compare an existing MTZ instead")
    ap.add_argument("--extra", default="", help="CONTROL lines before PROCESS")
    ap.add_argument("--keep", type=Path, help="work directory to keep")
    ap.add_argument("--json", type=Path, help="write metrics as JSON")
    args = ap.parse_args()

    if args.mtz:
        mtz = args.mtz
    else:
        work = args.keep or Path(tempfile.mkdtemp(prefix="bench_thaum_"))
        work.mkdir(parents=True, exist_ok=True)
        extra = args.extra.replace("\\n", "\n")
        if extra and not extra.endswith("\n"):
            extra += "\n"
        mtz = run_processinglab(
            args.exe, args.image, args.frames, args.nproc, work, extra
        )
        print(f"  work dir: {work}")
    report = compare(load(mtz), load(args.dials))
    print_report(report)
    if args.json:
        args.json.write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
