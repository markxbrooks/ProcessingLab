#!/usr/bin/env python3
"""Head-to-head ProcessingLab vs Kabsch XDS on one diffraction sweep.

Runs both pipelines on the same images / geometry seed, then extracts
indexing, integration and (optional) scaling metrics into a JSON summary
and a printed table.

Usage (from ProcessingLab root)::

  # Quick smoke test (first 50 frames)
  python3 processing-lab/tools/compare_xds_vs_processinglab.py \\
      --image ~/projects/76c1b_0002212/76c1b_0002212_1_0001.cbf \\
      --range 1:50 --workdir /tmp/xds_pl_cmp --aimless

  # Full 900-frame CBF sweep (matches the Oct 2026 canvas run)
  python3 processing-lab/tools/compare_xds_vs_processinglab.py \\
      --image ~/projects/76c1b_0002212/76c1b_0002212_1_0001.cbf \\
      --xds-inp processing-lab/tools/data/76c1b_xds_seed.INP \\
      --workdir ~/projects/76c1b_0002212/xds_vs_pl_rerun --aimless

  # Re-parse an existing workdir without re-running
  python3 processing-lab/tools/compare_xds_vs_processinglab.py \\
      --workdir ~/projects/76c1b_0002212/xds_vs_pl_20261025 --parse-only
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, MutableMapping, Optional, Sequence

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PL = ROOT / "bin" / "ProcessingLab"


def _first_existing(*candidates: Path | str | None) -> Path:
    """The first candidate that is a file (the last one if none is)."""
    paths = [Path(c) for c in candidates if c]
    return next((c for c in paths if c.is_file()), paths[-1])


# PATH first, then the usual macOS (Apple-silicon XDS, /Applications CCP4)
# and Linux locations.
DEFAULT_XDS = _first_existing(
    shutil.which("xds"),
    Path.home() / "Downloads" / "XDS-Apple_M1" / "xds",
    Path.home() / "Downloads" / "XDS-gfortran_Linux_x86_64" / "xds",
)
DEFAULT_IMAGE = (
    Path.home() / "projects" / "76c1b_0002212" / "76c1b_0002212_1_0001.cbf"
)
DEFAULT_POINTLESS = _first_existing(
    shutil.which("pointless"),
    Path(os.environ.get("CBIN", "/nonexistent")) / "pointless",
    Path("/Applications/ccp4-9/bin/pointless"),
    Path("/apps/ccp4/ccp4-9/bin/pointless"),
)
DEFAULT_AIMLESS = _first_existing(
    shutil.which("aimless"),
    Path(os.environ.get("CBIN", "/nonexistent")) / "aimless",
    Path("/Applications/ccp4-9/bin/aimless"),
    Path("/apps/ccp4/ccp4-9/bin/aimless"),
)

# Pilatus 6M geometry for 76c1b_0002212 (CBF header). Used when no --xds-inp.
DEFAULT_GEOMETRY = {
    "DETECTOR": "PILATUS",
    "MINIMUM_VALID_PIXEL_VALUE": "0",
    "OVERLOAD": "244849",
    "SENSOR_THICKNESS": "0.320",
    "NX": "2463",
    "NY": "2527",
    "QX": "0.172",
    "QY": "0.172",
    "ORGX": "1245.82",
    "ORGY": "1209.07",
    "DETECTOR_DISTANCE": "385.19",
    "OSCILLATION_RANGE": "0.2",
    "X-RAY_WAVELENGTH": "0.97949",
    "STARTING_ANGLE": "0.0",
    "STARTING_FRAME": "1",
    "ROTATION_AXIS": "1.0 0.0 0.0",
    "DIRECTION_OF_DETECTOR_X-AXIS": "1.0 0.0 0.0",
    "DIRECTION_OF_DETECTOR_Y-AXIS": "0.0 1.0 0.0",
    "INCIDENT_BEAM_DIRECTION": "0.0 0.0 1.0",
    "FRACTION_OF_POLARIZATION": "0.99",
    "POLARIZATION_PLANE_NORMAL": "0.0 1.0 0.0",
}

ICE_RINGS = (
    (3.93, 3.87),
    (3.70, 3.64),
    (3.47, 3.41),
    (2.70, 2.64),
    (2.28, 2.22),
    (2.102, 2.042),
    (1.94, 1.92),
    (1.88, 1.82),
    (1.74, 1.70),
)


def _now_stamp() -> str:
    """Return a compact local timestamp for default workdirs."""
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _run(
    cmd: Sequence[str],
    *,
    cwd: Path,
    log: Path,
    env: Optional[Mapping[str, str]] = None,
    stdin_text: str = "",
) -> float:
    """Run a command, tee stdout/stderr to *log*, return wall seconds.

    :param cmd: argv
    :param cwd: working directory
    :param log: log file path
    :param env: optional environment overlay
    :param stdin_text: keywords for the program's standard input (CCP4
        programs read keywords until end of file; empty = none)
    :return: elapsed wall-clock seconds
    :raises subprocess.CalledProcessError: on non-zero exit
    """
    log.parent.mkdir(parents=True, exist_ok=True)
    merged = os.environ.copy()
    if env:
        merged.update(env)
    print(f"  $ {' '.join(shlex.quote(c) for c in cmd)}", flush=True)
    t0 = time.perf_counter()
    with log.open("w") as fh:
        fh.write(f"# cwd={cwd}\n# cmd={' '.join(cmd)}\n\n")
        fh.flush()
        proc = subprocess.run(
            list(cmd),
            cwd=str(cwd),
            input=stdin_text,
            text=True,
            stdout=fh,
            stderr=subprocess.STDOUT,
            env=merged,
            check=False,
        )
    elapsed = time.perf_counter() - t0
    with log.open("a") as fh:
        fh.write(f"\n# exit={proc.returncode} elapsed={elapsed:.2f}s\n")
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)
    return elapsed


def _parse_range(text: str) -> tuple[int, int]:
    """Parse ``FIRST:LAST`` or ``FIRST-LAST`` into inclusive ints."""
    m = re.fullmatch(r"\s*(\d+)\s*[:\-]\s*(\d+)\s*", text)
    if not m:
        raise argparse.ArgumentTypeError(f"bad range {text!r}; want FIRST:LAST")
    first, last = int(m.group(1)), int(m.group(2))
    if last < first:
        raise argparse.ArgumentTypeError(f"range end < start: {text}")
    return first, last


def _cbf_template(image: Path) -> str:
    """Build an XDS ``NAME_TEMPLATE_OF_DATA_FRAMES`` for a numbered CBF.

    :param image: one frame path (``…_0001.cbf``)
    :return: absolute template with ``????`` digits
    """
    name = image.name
    m = re.search(r"(\d+)(\.cbf)$", name, re.IGNORECASE)
    if not m:
        raise ValueError(f"cannot derive CBF template from {image}")
    digits = m.group(1)
    stem = name[: m.start(1)] + ("?" * len(digits)) + m.group(2)
    return str((image.parent / stem).resolve())


def write_xds_inp(
    path: Path,
    *,
    image: Path,
    data_range: tuple[int, int],
    spot_range: tuple[int, int],
    job: str,
    nproc: int,
    space_group: int,
    seed: Optional[Path],
    include_resolution: tuple[float, float],
    exclude_ice: bool,
) -> None:
    """Write an ``XDS.INP`` for the comparison run.

    :param path: output path
    :param image: first-frame / master path used to build the template
    :param data_range: inclusive frame range
    :param spot_range: COLSPOT / FINDSPOTS range
    :param job: XDS ``JOB=`` value
    :param nproc: ``MAXIMUM_NUMBER_OF_PROCESSORS``
    :param space_group: ``SPACE_GROUP_NUMBER`` (0 = unknown)
    :param seed: optional existing XDS.INP to copy geometry from
    :param include_resolution: ``INCLUDE_RESOLUTION_RANGE``
    :param exclude_ice: add standard ice-ring excludes
    """
    lines: list[str] = [
        f"! compare_xds_vs_processinglab.py  {datetime.now(timezone.utc).isoformat()}",
        f"JOB= {job}",
    ]
    if seed is not None:
        for raw in seed.read_text().splitlines():
            s = raw.strip()
            if not s or s.startswith("!"):
                continue
            key = s.split("=", 1)[0].strip().upper()
            if key in {
                "JOB",
                "NAME_TEMPLATE_OF_DATA_FRAMES",
                "DATA_RANGE",
                "SPOT_RANGE",
                "BACKGROUND_RANGE",
                "MAXIMUM_NUMBER_OF_PROCESSORS",
                "MAXIMUM_NUMBER_OF_PROCESSOR",
                "SPACE_GROUP_NUMBER",
                "INCLUDE_RESOLUTION_RANGE",
                "STRONG_PIXEL",  # obsolete in recent XDS
            }:
                continue
            if key.startswith("EXCLUDE_RESOLUTION_RANGE") and exclude_ice:
                continue
            lines.append(raw.rstrip())
    else:
        g = DEFAULT_GEOMETRY
        lines.extend(
            [
                f"DETECTOR={g['DETECTOR']} "
                f"MINIMUM_VALID_PIXEL_VALUE={g['MINIMUM_VALID_PIXEL_VALUE']} "
                f"OVERLOAD={g['OVERLOAD']}",
                f"SENSOR_THICKNESS={g['SENSOR_THICKNESS']}",
                f"NX={g['NX']} NY={g['NY']} QX={g['QX']} QY={g['QY']}",
                f"ORGX={g['ORGX']} ORGY={g['ORGY']}",
                f"DETECTOR_DISTANCE={g['DETECTOR_DISTANCE']}",
                f"OSCILLATION_RANGE={g['OSCILLATION_RANGE']}",
                f"X-RAY_WAVELENGTH={g['X-RAY_WAVELENGTH']}",
                f"STARTING_ANGLE={g['STARTING_ANGLE']} "
                f"STARTING_FRAME={g['STARTING_FRAME']}",
                f"ROTATION_AXIS={g['ROTATION_AXIS']}",
                f"DIRECTION_OF_DETECTOR_X-AXIS={g['DIRECTION_OF_DETECTOR_X-AXIS']}",
                f"DIRECTION_OF_DETECTOR_Y-AXIS={g['DIRECTION_OF_DETECTOR_Y-AXIS']}",
                f"INCIDENT_BEAM_DIRECTION={g['INCIDENT_BEAM_DIRECTION']}",
                f"FRACTION_OF_POLARIZATION={g['FRACTION_OF_POLARIZATION']}",
                f"POLARIZATION_PLANE_NORMAL={g['POLARIZATION_PLANE_NORMAL']}",
                "TRUSTED_REGION=0.0 1.41",
                "MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT=3",
                "SEPMIN=6.0",
                "CLUSTER_RADIUS=2",
                "VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS=7000 30000",
                "UNIT_CELL_CONSTANTS=0 0 0 0 0 0",
                "REFINE(IDXREF)=BEAM AXIS ORIENTATION CELL",
                "REFINE(INTEGRATE)=POSITION BEAM ORIENTATION CELL",
                "REFINE(CORRECT)=DISTANCE BEAM AXIS ORIENTATION CELL",
            ]
        )

    template = _cbf_template(image) if image.suffix.lower() == ".cbf" else str(image.resolve())
    first, last = data_range
    s0, s1 = spot_range
    bg1 = min(10, last)
    lines.extend(
        [
            f"NAME_TEMPLATE_OF_DATA_FRAMES={template}",
            f"DATA_RANGE={first} {last}",
            f"SPOT_RANGE={s0} {s1}",
            f"BACKGROUND_RANGE={first} {bg1}",
            f"MAXIMUM_NUMBER_OF_PROCESSORS={nproc}",
            f"SPACE_GROUP_NUMBER={space_group}",
            f"INCLUDE_RESOLUTION_RANGE={include_resolution[0]} {include_resolution[1]}",
        ]
    )
    if exclude_ice:
        for lo, hi in ICE_RINGS:
            lines.append(f"EXCLUDE_RESOLUTION_RANGE={lo} {hi}")

    path.write_text("\n".join(lines) + "\n")


def run_processinglab(
    *,
    exe: Path,
    image: Path,
    workdir: Path,
    data_range: tuple[int, int],
    xds_inp: Path,
    nproc: int,
    refine: str,
    exclude_ice: bool,
) -> float:
    """Run ``ProcessingLab run`` into *workdir*.

    :return: wall-clock seconds
    """
    if workdir.exists():
        raise FileExistsError(f"PL output dir already exists: {workdir}")
    first, last = data_range
    cmd = [
        str(exe),
        "run",
        "--image",
        str(image.resolve()),
        "--range",
        f"{first}:{last}",
        "--xds-inp",
        str(xds_inp.resolve()),
        "--refine",
        refine,
        "--output-dir",
        str(workdir.resolve()),
    ]
    if exclude_ice:
        cmd.append("--exclude-ice")
    env = {"OMP_NUM_THREADS": str(nproc)}
    # ProcessingLab writes its own process.log inside workdir; keep a
    # wrapper log next to it for the launcher argv / timing.
    return _run(cmd, cwd=workdir.parent, log=workdir.parent / "pl_launch.log", env=env)


def run_xds(
    *,
    exe: Path,
    xds_dir: Path,
    image: Path,
    data_range: tuple[int, int],
    spot_range: tuple[int, int],
    seed: Optional[Path],
    nproc: int,
    space_group_index: int,
    space_group_correct: int,
    include_resolution: tuple[float, float],
    exclude_ice: bool,
    two_pass: bool,
) -> dict[str, float]:
    """Run XDS indexing then integration/correct.

    :param space_group_index: SG for IDXREF (usually 0)
    :param space_group_correct: SG for CORRECT (0 keeps IDXREF choice / P1)
    :param two_pass: if True, IDXREF then INTEGRATE+CORRECT separately
    :return: mapping of stage name → wall seconds
    """
    xds_dir.mkdir(parents=True, exist_ok=True)
    times: dict[str, float] = {}

    def write(job: str, sg: int) -> Path:
        inp = xds_dir / "XDS.INP"
        write_xds_inp(
            inp,
            image=image,
            data_range=data_range,
            spot_range=spot_range,
            job=job,
            nproc=nproc,
            space_group=sg,
            seed=seed,
            include_resolution=include_resolution,
            exclude_ice=exclude_ice,
        )
        return inp

    if two_pass:
        write("XYCORR INIT COLSPOT IDXREF", space_group_index)
        times["idxref"] = _run(
            [str(exe)], cwd=xds_dir, log=xds_dir / "xds_idxref.log"
        )
        # Keep refined XPARM; switch JOB for integrate.
        write("DEFPIX INTEGRATE CORRECT", space_group_correct)
        times["integrate_correct"] = _run(
            [str(exe)], cwd=xds_dir, log=xds_dir / "xds_integrate.log"
        )
    else:
        write(
            "XYCORR INIT COLSPOT IDXREF DEFPIX INTEGRATE CORRECT",
            space_group_index if space_group_correct == 0 else space_group_correct,
        )
        # If index SG and correct SG differ, force two-pass.
        if space_group_correct not in (0, space_group_index):
            write("XYCORR INIT COLSPOT IDXREF", space_group_index)
            times["idxref"] = _run(
                [str(exe)], cwd=xds_dir, log=xds_dir / "xds_idxref.log"
            )
            write("DEFPIX INTEGRATE CORRECT", space_group_correct)
            times["integrate_correct"] = _run(
                [str(exe)], cwd=xds_dir, log=xds_dir / "xds_integrate.log"
            )
        else:
            times["all"] = _run(
                [str(exe)], cwd=xds_dir, log=xds_dir / "xds_all.log"
            )
    return times


def run_aimless(
    *,
    hkl: Path,
    out_dir: Path,
    pointless: Path,
    aimless: Path,
    xds_ascii: bool = False,
    space_group: Optional[str] = None,
    resolution: Optional[float] = None,
) -> float:
    """Run Pointless then Aimless on an unmerged MTZ or XDS_ASCII.HKL.

    :param hkl: Unmerged MTZ (``HKLIN``) or XDS ``XDS_ASCII.HKL`` (``XDSIN``).
    :param xds_ascii: *hkl* is an XDS_ASCII file.
    :param space_group: Force this space group in Pointless (e.g. ``"P 1 2 1"``),
        so both programs' data are scaled in the same group.
    :param resolution: Aimless high-resolution cutoff (Å), or ``None``.
    :return: wall-clock seconds
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    pointless_mtz = out_dir / "pointless.mtz"
    aimless_mtz = out_dir / "aimless.mtz"
    t0 = time.perf_counter()
    # CCP4 logical names are separate arguments (``HKLIN file``, not
    # ``HKLIN=file``), and the programs read keywords from stdin until EOF.
    _run(
        [str(pointless), "XDSIN" if xds_ascii else "HKLIN", str(hkl), "HKLOUT", str(pointless_mtz)],
        cwd=out_dir,
        log=out_dir / "pointless.log",
        stdin_text=f"SPACEGROUP {space_group}\n" if space_group else "",
    )
    _run(
        [str(aimless), "HKLIN", str(pointless_mtz), "HKLOUT", str(aimless_mtz)],
        cwd=out_dir,
        log=out_dir / "aimless.log",
        stdin_text=f"RESOLUTION HIGH {resolution}\n" if resolution else "",
    )
    return time.perf_counter() - t0


def _search(text: str, pattern: str, flags: int = 0) -> Optional[re.Match[str]]:
    return re.search(pattern, text, flags)


def parse_pl(pl_dir: Path) -> dict[str, Any]:
    """Extract metrics from a ProcessingLab run directory."""
    out: dict[str, Any] = {"dir": str(pl_dir)}
    md = pl_dir / "process.md"
    if md.is_file():
        text = md.read_text()
        out["process_md"] = str(md)
        m = _search(text, r"Finished:.*?\((\d+(?:\.\d+)?) s\)")
        if m:
            out["wall_s"] = float(m.group(1))
        m = _search(
            text,
            r"Lattice:\s*(\S+),\s*space group\s*(\d+)",
        )
        if m:
            out["lattice"] = m.group(1)
            out["space_group"] = int(m.group(2))
        m = _search(
            text,
            r"Cell:\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*Å,\s*"
            r"([\d.]+)\s+([\d.]+)\s+([\d.]+)°",
        )
        if m:
            out["cell"] = [float(x) for x in m.groups()]
        # Indexing block (earlier in md) vs refined geometry
        cells = re.findall(
            r"Cell:\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*Å,\s*"
            r"([\d.]+)\s+([\d.]+)\s+([\d.]+)°",
            text,
        )
        if cells:
            out["cell_index"] = [float(x) for x in cells[0]]
            out["cell_refined"] = [float(x) for x in cells[-1]]
        m = _search(text, r"Indexed\s+(\d+)%\s+of\s+(\d+)\s+spots")
        if m:
            out["indexed_pct"] = int(m.group(1))
            out["nspots"] = int(m.group(2))
        m = _search(text, r"Mosaicity:\s*([\d.]+)°")
        if m:
            out["mosaicity"] = float(m.group(1))
        m = _search(text, r"integrated=(\d+)")
        if m:
            out["n_integrated"] = int(m.group(1))
        m = _search(
            text,
            r"\|\s*\*\*total\*\*\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*"
            r"([\d.]+)%\s*\|\s*([\d.]+)%\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)",
        )
        if m:
            out["unscaled"] = {
                "observed": int(m.group(1)),
                "unique": int(m.group(2)),
                "completeness_pct": float(m.group(3)),
                "rmeas_pct": float(m.group(4)),
                "isigi": float(m.group(5)),
                "cc12": float(m.group(6)),
            }
    pj = pl_dir / "process.json"
    if pj.is_file():
        try:
            out["process_json"] = json.loads(pj.read_text())
        except json.JSONDecodeError:
            pass
    mtz = pl_dir / "process.mtz"
    if mtz.is_file():
        out["mtz"] = str(mtz)
    return out


def parse_xds(xds_dir: Path) -> dict[str, Any]:
    """Extract metrics from an XDS working directory."""
    out: dict[str, Any] = {"dir": str(xds_dir)}

    idx = xds_dir / "IDXREF.LP"
    if idx.is_file():
        text = idx.read_text()
        m = _search(
            text,
            r"(\d+)\s+OUT OF\s+(\d+)\s+SPOTS INDEXED",
        )
        if m:
            indexed, total = int(m.group(1)), int(m.group(2))
            out["indexed"] = indexed
            out["spots"] = total
            out["indexed_pct"] = round(100.0 * indexed / total, 1) if total else None
        # Last refined cell in IDXREF
        cells = re.findall(
            r"UNIT CELL PARAMETERS\s+"
            r"([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+"
            r"([\d.]+)\s+([\d.]+)\s+([\d.]+)",
            text,
        )
        if cells:
            out["cell_index"] = [float(x) for x in cells[-1]]
        m = _search(text, r"CRYSTAL MOSAICITY \(DEGREES\)\s+([\d.]+)")
        if m:
            out["mosaicity_index"] = float(m.group(1))

    gx = xds_dir / "GXPARM.XDS"
    if gx.is_file():
        # Line 4 of Kabsch GXPARM: sg a b c al be ga
        lines = [ln for ln in gx.read_text().splitlines() if ln.strip()]
        if len(lines) >= 4:
            parts = lines[3].split()
            if len(parts) >= 7:
                out["space_group"] = int(float(parts[0]))
                out["cell_refined"] = [float(x) for x in parts[1:7]]

    integ = xds_dir / "INTEGRATE.LP"
    if integ.is_file():
        text = integ.read_text()
        mosaics = re.findall(
            r"CRYSTAL MOSAICITY \(DEGREES\)\s+([\d.]+)", text
        )
        if mosaics:
            out["mosaicity"] = float(mosaics[-1])
        m = _search(text, r"(\d+)\s+REFLECTIONS SAVED ON FILE")
        if m:
            out["n_integrated"] = int(m.group(1))

    corr = xds_dir / "CORRECT.LP"
    if corr.is_file():
        text = corr.read_text()
        # Last "total" row in a COMPLETENESS AND QUALITY table
        totals = re.findall(
            r"^\s*total\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)%\s+"
            r"([\d.]+)%\s+([\d.]+)%\s+(\d+)\s+([\d.]+)\s+"
            r"([\d.]+)%\s+([\d.]+)",
            text,
            flags=re.MULTILINE,
        )
        if totals:
            t = totals[-1]
            out["correct"] = {
                "observed": int(t[0]),
                "unique": int(t[1]),
                "possible": int(t[2]),
                "completeness_pct": float(t[3]),
                "rfactor_obs_pct": float(t[4]),
                "rfactor_exp_pct": float(t[5]),
                "compared": int(t[6]),
                "isigi": float(t[7]),
                "rmeas_pct": float(t[8]),
                "cc12": float(t[9]),
            }
        m = _search(text, r"ISa\s*\n\s*[\d.E+-]+\s+[\d.E+-]+\s+([\d.]+)")
        if not m:
            m = _search(
                text,
                r"a\s+b\s+ISa\s*\n\s*[\d.E+-]+\s+[\d.E+-]+\s+([\d.]+)",
            )
        if m:
            out["isa"] = float(m.group(1))
        m = _search(
            text,
            r"SELECTED SPACE GROUP AND UNIT CELL.*?"
            r"SPACE_GROUP_NUMBER=\s*(\d+)\s*"
            r"UNIT_CELL_CONSTANTS=\s*"
            r"([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+"
            r"([\d.]+)\s+([\d.]+)\s+([\d.]+)",
            flags=re.DOTALL,
        )
        if m:
            out["correct_space_group"] = int(m.group(1))
            out["correct_cell"] = [float(x) for x in m.groups()[1:]]

    wall_parts: dict[str, float] = {}
    for name in ("xds_all.log", "xds_integrate.log", "xds_idxref.log"):
        log = xds_dir / name
        if not log.is_file():
            continue
        text = log.read_text()
        # Prefer GNU time / wrapper line written by this script.
        m = _search(text, r"^# exit=\d+ elapsed=([\d.]+)s", flags=re.MULTILINE)
        if not m:
            m = _search(
                text,
                r"^elapsed ([\d.]+) s\s*$",
                flags=re.MULTILINE,
            )
        if not m:
            m = _search(
                text,
                r"Total elapsed wall-clock time for XDS\s+([\d.]+) sec",
            )
        if m:
            wall_parts[name] = float(m.group(1))
    if wall_parts:
        out["wall_s_parts"] = wall_parts
        out["wall_s"] = sum(wall_parts.values())

    return out


def _first_log(directory: Path, *patterns: str) -> Optional[Path]:
    """Return the first existing log matching *patterns* (name or glob)."""
    for pat in patterns:
        direct = directory / pat
        if direct.is_file():
            return direct
        matches = sorted(directory.glob(pat))
        if matches:
            return matches[0]
    return None


def parse_aimless(ccp4_dir: Path) -> dict[str, Any]:
    """Extract overall Aimless / Pointless metrics."""
    out: dict[str, Any] = {"dir": str(ccp4_dir)}
    plog = _first_log(ccp4_dir, "pointless.log", "*pointless*.log")
    if plog is not None:
        text = plog.read_text()
        m = _search(text, r"Best Solution:\s*space group\s+(.+)")
        if m:
            out["space_group"] = m.group(1).strip()
    alog = _first_log(ccp4_dir, "aimless.log", "*aimless*.log")
    if alog is None:
        return out
    text = alog.read_text()
    out["log"] = str(alog)

    # Analysis against resolution: Rmrg Rfull Rcum Rmeas Rpim Nmeas ...
    m = _search(
        text,
        r"^Overall:\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+"
        r"([\d.]+)\s+(\d+)\s+",
        flags=re.MULTILINE,
    )
    if m:
        out["rmeas"] = float(m.group(4))
        out["rpim"] = float(m.group(5))
        out["nmeas"] = int(m.group(6))

    # Completeness table: Nmeas Nref Ncent %poss C%poss Mlplct ...
    m = _search(
        text,
        r"^Overall:\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
        flags=re.MULTILINE,
    )
    if m:
        out["observed"] = int(m.group(1))
        out["unique"] = int(m.group(2))
        out["completeness_pct"] = float(m.group(4))
        out["multiplicity"] = float(m.group(6))

    # CC(1/2) table header mentions CC1/2; Overall: CCanom Nanom RCRanom CC1/2 ...
    m = _search(
        text,
        r"\$TABLE:\s*Correlations CC\(1/2\).*?^Overall:\s+"
        r"[-\d.]+\s+\d+\s+[\d.]+\s+([\d.]+)",
        flags=re.DOTALL | re.MULTILINE,
    )
    if m:
        out["cc12"] = float(m.group(1))

    # Mn(I/sd) on the Rmeas-vs-resolution Overall line:
    # Rmrg Rfull Rcum Rmeas Rpim Nmeas AvI RMSdev sd I/RMS Mn(I/sd) ...
    if "rmeas" in out and "nmeas" in out:
        for line in text.splitlines():
            if not line.startswith("Overall:"):
                continue
            parts = line.split()
            if len(parts) < 12:
                continue
            try:
                if int(parts[6]) == out["nmeas"] and abs(float(parts[4]) - out["rmeas"]) < 1e-6:
                    out["isigi"] = float(parts[11])
                    break
            except ValueError:
                continue
    return out


def _fmt_cell(cell: Optional[Sequence[float]]) -> str:
    if not cell or len(cell) < 6:
        return "—"
    a, b, c, al, be, ga = cell[:6]
    return f"{a:.2f} {b:.2f} {c:.2f} / {al:.2f} {be:.2f} {ga:.2f}"


def print_summary(summary: Mapping[str, Any]) -> None:
    """Print a compact side-by-side report."""
    pl = summary.get("processinglab") or {}
    xds = summary.get("xds") or {}
    aim = summary.get("aimless") or {}
    aim_x = summary.get("aimless_xds") or {}

    rows = [
        ("Index cell", _fmt_cell(pl.get("cell_index")), _fmt_cell(xds.get("cell_index"))),
        ("Refined cell", _fmt_cell(pl.get("cell_refined") or pl.get("cell")),
         _fmt_cell(xds.get("cell_refined"))),
        ("% indexed", pl.get("indexed_pct", "—"), xds.get("indexed_pct", "—")),
        ("Mosaicity °", pl.get("mosaicity", "—"), xds.get("mosaicity", "—")),
        ("N integrated", pl.get("n_integrated", "—"), xds.get("n_integrated", "—")),
        ("Wall s", pl.get("wall_s", "—"), xds.get("wall_s", "—")),
    ]
    xu = (xds.get("correct") or {})
    pu = (pl.get("unscaled") or {})
    rows.extend(
        [
            ("Rmeas %", pu.get("rmeas_pct", "—"), xu.get("rmeas_pct", "—")),
            ("CC½", pu.get("cc12", "—"), xu.get("cc12", "—")),
            ("I/σ", pu.get("isigi", "—"), xu.get("isigi", "—")),
            ("Completeness %", pu.get("completeness_pct", "—"),
             xu.get("completeness_pct", "—")),
        ]
    )

    print("\n=== ProcessingLab vs XDS ===")
    print(f"{'metric':<18} {'ProcessingLab':<36} {'XDS':<36}")
    print("-" * 90)
    for name, a, b in rows:
        print(f"{name:<18} {str(a):<36} {str(b):<36}")

    if aim or aim_x:
        print("\nPointless + Aimless (the same programs and settings for both):")
        print(f"{'metric':<18} {'ProcessingLab':<36} {'XDS (XDS_ASCII.HKL)':<36}")
        print("-" * 90)
        for k in ("space_group", "rmeas", "cc12", "isigi", "completeness_pct",
                  "multiplicity", "observed", "unique"):
            if k in aim or k in aim_x:
                print(f"{k:<18} {str(aim.get(k, '—')):<36} {str(aim_x.get(k, '—')):<36}")
    print()


def build_argparser() -> argparse.ArgumentParser:
    """Create the CLI parser."""
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\n".join(__doc__.splitlines()[1:]),
    )
    ap.add_argument(
        "--image",
        type=Path,
        default=DEFAULT_IMAGE,
        help="first CBF frame or HDF5 master (default: 76c1b_0001.cbf)",
    )
    ap.add_argument(
        "--range",
        dest="frame_range",
        type=_parse_range,
        default=None,
        help="inclusive frame range FIRST:LAST (default: 1:900)",
    )
    ap.add_argument(
        "--spot-range",
        type=_parse_range,
        default=None,
        help="spot-finding range (default: 1:min(50, last))",
    )
    ap.add_argument(
        "--workdir",
        type=Path,
        default=None,
        help="comparison root (default: ./xds_vs_pl_<timestamp>)",
    )
    ap.add_argument("--pl-exe", type=Path, default=DEFAULT_PL)
    ap.add_argument("--xds-exe", type=Path, default=DEFAULT_XDS)
    ap.add_argument(
        "--xds-inp",
        type=Path,
        default=None,
        help="geometry seed XDS.INP (default: built-in 76c1b Pilatus params)",
    )
    ap.add_argument("--nproc", type=int, default=4)
    ap.add_argument(
        "--refine",
        default="orientation,cell,distance,beam,tilt_twist,mosaicity",
        help="ProcessingLab --refine list",
    )
    ap.add_argument(
        "--dmin",
        type=float,
        default=1.5,
        help="high-resolution limit for XDS INCLUDE_RESOLUTION_RANGE",
    )
    ap.add_argument(
        "--dmax",
        type=float,
        default=50.0,
        help="low-resolution limit for XDS INCLUDE_RESOLUTION_RANGE",
    )
    ap.add_argument(
        "--space-group-index",
        type=int,
        default=0,
        help="XDS SPACE_GROUP_NUMBER for IDXREF (0=unknown)",
    )
    ap.add_argument(
        "--space-group-correct",
        type=int,
        default=0,
        help="XDS SPACE_GROUP_NUMBER for CORRECT (0=leave / P1 path)",
    )
    ap.add_argument(
        "--exclude-ice",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="exclude standard ice rings (default: on)",
    )
    ap.add_argument(
        "--two-pass-xds",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="run IDXREF then INTEGRATE+CORRECT separately (default: on)",
    )
    ap.add_argument("--skip-pl", action="store_true")
    ap.add_argument("--skip-xds", action="store_true")
    ap.add_argument(
        "--aimless",
        action="store_true",
        help="run Pointless+Aimless on ProcessingLab's process.mtz and XDS's XDS_ASCII.HKL",
    )
    ap.add_argument(
        "--aimless-spacegroup",
        default=None,
        help="force this space group in Pointless for both (e.g. 'P 1 2 1')",
    )
    ap.add_argument(
        "--aimless-resolution",
        type=float,
        default=None,
        help="Aimless high-resolution cutoff (Å) for both",
    )
    ap.add_argument("--pointless", type=Path, default=DEFAULT_POINTLESS)
    ap.add_argument("--aimless-exe", type=Path, default=DEFAULT_AIMLESS)
    ap.add_argument(
        "--parse-only",
        action="store_true",
        help="only parse --workdir (expects pl/ and xds/ subdirs, or a bare XDS dir)",
    )
    ap.add_argument("--json", type=Path, default=None, help="summary JSON path")
    return ap


def resolve_dirs(workdir: Path, parse_only: bool) -> tuple[Path, Path, Path]:
    """Return ``(workdir, pl_dir, xds_dir)``, creating layout as needed."""
    pl_dir = workdir / "pl"
    xds_dir = workdir / "xds"
    if parse_only:
        # Accept a bare prior XDS workdir (no pl/ child).
        if not pl_dir.is_dir() and (workdir / "IDXREF.LP").is_file():
            xds_dir = workdir
        if not pl_dir.is_dir():
            # Try sibling processinglab session
            pass
        return workdir, pl_dir, xds_dir
    workdir.mkdir(parents=True, exist_ok=True)
    return workdir, pl_dir, xds_dir


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entry point.

    :param argv: optional argv override
    :return: process exit code
    """
    ap = build_argparser()
    args = ap.parse_args(argv)

    workdir = args.workdir or Path(f"xds_vs_pl_{_now_stamp()}").resolve()
    workdir, pl_dir, xds_dir = resolve_dirs(workdir, args.parse_only)
    data_range = args.frame_range or (1, 900)
    spot_range = args.spot_range or (1, min(50, data_range[1]))
    include_resolution = (args.dmax, args.dmin)

    summary: MutableMapping[str, Any] = {
        "created": datetime.now(timezone.utc).isoformat(),
        "image": str(args.image.resolve()) if args.image.exists() else str(args.image),
        "range": list(data_range),
        "spot_range": list(spot_range),
        "workdir": str(workdir.resolve()),
    }

    seed_inp = workdir / "geometry_seed.XDS.INP"
    if not args.parse_only:
        if not args.image.is_file():
            print(f"error: image not found: {args.image}", file=sys.stderr)
            return 2
        write_xds_inp(
            seed_inp,
            image=args.image,
            data_range=data_range,
            spot_range=spot_range,
            job="XYCORR INIT COLSPOT IDXREF",
            nproc=args.nproc,
            space_group=args.space_group_index,
            seed=args.xds_inp,
            include_resolution=include_resolution,
            exclude_ice=args.exclude_ice,
        )
        summary["geometry_seed"] = str(seed_inp)

        if not args.skip_pl:
            if not args.pl_exe.is_file():
                print(f"error: ProcessingLab not found: {args.pl_exe}", file=sys.stderr)
                return 2
            print("== ProcessingLab ==")
            try:
                pl_s = run_processinglab(
                    exe=args.pl_exe,
                    image=args.image,
                    workdir=pl_dir,
                    data_range=data_range,
                    xds_inp=seed_inp,
                    nproc=args.nproc,
                    refine=args.refine,
                    exclude_ice=args.exclude_ice,
                )
                summary["pl_wall_s"] = pl_s
            except (subprocess.CalledProcessError, FileExistsError) as exc:
                print(f"error: ProcessingLab failed: {exc}", file=sys.stderr)
                return 1

        if not args.skip_xds:
            if not args.xds_exe.is_file():
                print(f"error: xds not found: {args.xds_exe}", file=sys.stderr)
                return 2
            print("== XDS ==")
            try:
                xds_times = run_xds(
                    exe=args.xds_exe,
                    xds_dir=xds_dir,
                    image=args.image,
                    data_range=data_range,
                    spot_range=spot_range,
                    seed=args.xds_inp or seed_inp,
                    nproc=args.nproc,
                    space_group_index=args.space_group_index,
                    space_group_correct=args.space_group_correct,
                    include_resolution=include_resolution,
                    exclude_ice=args.exclude_ice,
                    two_pass=args.two_pass_xds,
                )
                summary["xds_wall_s"] = xds_times
            except subprocess.CalledProcessError as exc:
                print(f"error: XDS failed: {exc}", file=sys.stderr)
                return 1

        if args.aimless and not args.skip_pl:
            mtz = pl_dir / "process.mtz"
            if not mtz.is_file():
                print(f"warning: no MTZ at {mtz}; skipping Aimless", file=sys.stderr)
            elif not args.pointless.is_file() or not args.aimless_exe.is_file():
                print("warning: pointless/aimless not found; skipping", file=sys.stderr)
            else:
                print("== Pointless + Aimless ==")
                ccp4_dir = workdir / "ccp4"
                try:
                    summary["aimless_wall_s"] = run_aimless(
                        hkl=mtz,
                        out_dir=ccp4_dir,
                        pointless=args.pointless,
                        aimless=args.aimless_exe,
                        space_group=args.aimless_spacegroup,
                        resolution=args.aimless_resolution,
                    )
                except subprocess.CalledProcessError as exc:
                    print(f"warning: Aimless failed: {exc}", file=sys.stderr)

        ascii_hkl = xds_dir / "XDS_ASCII.HKL"
        if args.aimless and not args.skip_xds and ascii_hkl.is_file() and args.aimless_exe.is_file():
            print("== Pointless + Aimless (XDS) ==")
            try:
                summary["aimless_xds_wall_s"] = run_aimless(
                    hkl=ascii_hkl,
                    out_dir=workdir / "ccp4_xds",
                    pointless=args.pointless,
                    aimless=args.aimless_exe,
                    xds_ascii=True,
                    space_group=args.aimless_spacegroup,
                    resolution=args.aimless_resolution,
                )
            except subprocess.CalledProcessError as exc:
                print(f"warning: Aimless (XDS) failed: {exc}", file=sys.stderr)

    # Parse whatever exists
    if pl_dir.is_dir():
        summary["processinglab"] = parse_pl(pl_dir)
        # Split-session layout: index metrics live in sibling 02-index
        # (resolve symlinks so pl -> …/03-integrate sees …/02-index).
        pl_real = pl_dir.resolve()
        for sibling in (
            pl_real.parent / "02-index",
            pl_dir.parent / "02-index",
            workdir / "02-index",
        ):
            if sibling.is_dir() and sibling.resolve() != pl_real:
                idx = parse_pl(sibling)
                pl_sum = summary["processinglab"]
                if idx.get("cell"):
                    pl_sum["cell_index"] = idx["cell"]
                for key in ("indexed_pct", "nspots", "lattice", "space_group"):
                    if key in idx:
                        pl_sum.setdefault(key, idx[key])
                break
    elif args.workdir and (Path(args.workdir) / "03-integrate").is_dir():
        session = Path(args.workdir)
        summary["processinglab"] = parse_pl(session / "03-integrate")
        if (session / "02-index").is_dir():
            idx = parse_pl(session / "02-index")
            pl_sum = summary["processinglab"]
            if idx.get("cell"):
                pl_sum["cell_index"] = idx["cell"]
            for key in ("indexed_pct", "nspots", "lattice", "space_group"):
                if key in idx:
                    pl_sum.setdefault(key, idx[key])

    if xds_dir.is_dir():
        summary["xds"] = parse_xds(xds_dir)

    ccp4_dir = workdir / "ccp4"
    if not _first_log(ccp4_dir, "aimless.log", "*aimless*.log"):
        pl_real = pl_dir.resolve() if pl_dir.exists() else pl_dir
        for alt in (workdir / "04-ccp4", pl_real.parent / "04-ccp4"):
            if _first_log(alt, "aimless.log", "*aimless*.log"):
                ccp4_dir = alt
                break
    if ccp4_dir.is_dir():
        summary["aimless"] = parse_aimless(ccp4_dir)
    if (workdir / "ccp4_xds").is_dir():
        summary["aimless_xds"] = parse_aimless(workdir / "ccp4_xds")

    print_summary(summary)

    json_path = args.json or (workdir / "summary.json")
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(summary, indent=2, default=str))
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
