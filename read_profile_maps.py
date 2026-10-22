#!/usr/bin/env python3
"""Plot the average 3D spot profile from a ProcessingLab log.

Reads a ``PROFILE_MAP`` block when the log has one. Otherwise it reads the
older three-across dump, where each value is three characters wide and a
peak of 100 is glued to its neighbor (``46100`` is 46 then 100). The last
profile in the file is the one that is plotted.
"""

import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

_LEGACY_HEADER = re.compile(
    r"images\s+(-?\d+)\s+\.\.\.\s+\+?(-?\d+)"
    r"\s+from the calculated phi;\s+(\d+)x(\d+)\s+pixels"
)


def _resolve_log(argv):
    """Find the log. ``Processing.Log`` is accepted as ``Process.Log``.

    A relative name is tried in the working directory and beside this script.
    """
    name = argv[1] if len(argv) > 1 else "Process.Log"
    path = Path(name)
    candidates = [path]
    if not path.is_absolute():
        candidates.append(Path(__file__).resolve().parent / path)
    if path.name.lower() == "processing.log":
        renamed = path.with_name("Process.Log")
        candidates.append(renamed)
        if not renamed.is_absolute():
            candidates.append(Path(__file__).resolve().parent / renamed.name)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return path


def _read_marked(lines):
    maps = []
    i = 0
    while i < len(lines):
        if not lines[i].startswith("PROFILE_MAP "):
            i += 1
            continue
        meta = dict(tok.split("=", 1) for tok in lines[i].split()[1:])
        nphi, side = int(meta["nphi"]), int(meta["side"])
        vol = np.zeros((nphi, side, side), dtype=int)
        i += 1
        ip = 0
        while i < len(lines) and lines[i] != "PROFILE_MAP_END":
            if lines[i].startswith("# phi="):
                rows = []
                i += 1
                while i < len(lines) and lines[i][:1].isspace() and lines[i].strip():
                    rows.append([int(x) for x in lines[i].split()])
                    i += 1
                vol[ip] = rows
                ip += 1
            else:
                i += 1
        maps.append(vol)
        i += 1
    return maps


def _fixed_panels(line, side, n_panels):
    """Split one legacy row into ``n_panels`` lists of ``side`` integers.

    Each panel is two spaces plus ``side`` fields of width 3.
    """
    panels = []
    i = 0
    for _ in range(n_panels):
        i += 2
        row = []
        for _col in range(side):
            row.append(int(line[i:i + 3]))
            i += 3
        panels.append(row)
    return panels


def _read_legacy(lines):
    maps = []
    i = 0
    while i < len(lines):
        match = _LEGACY_HEADER.search(lines[i])
        i += 1
        if match is None:
            continue
        phi_lo = int(match.group(1))
        phi_hi = int(match.group(2))
        side = int(match.group(3))
        nphi = phi_hi - phi_lo + 1
        vol = np.zeros((nphi, side, side), dtype=int)
        ip = 0
        while ip < nphi and i < len(lines):
            n_panels = min(3, nphi - ip)
            width = n_panels * (2 + side * 3)
            if len(lines[i]) != width:
                i += 1
                continue
            for _row in range(side):
                if i >= len(lines) or len(lines[i]) != width:
                    break
                for panel, row in enumerate(_fixed_panels(lines[i], side, n_panels)):
                    vol[ip + panel, _row] = row
                i += 1
            ip += n_panels
        maps.append(vol)
    return maps


def read_profile_maps(path):
    """Return every average profile in ``path``, oldest first.

    Each array has shape ``(nphi, side, side)``. Section 0 is the most
    negative image offset.

    :param path: Log file.
    :returns: One array per profile. Empty when the log has none.
    """
    lines = Path(path).read_text(errors="replace").splitlines()
    marked = _read_marked(lines)
    if marked:
        return marked
    return _read_legacy(lines)


def main():
    path = _resolve_log(sys.argv)
    if not path.is_file():
        sys.exit(f"no log at {path}")
    maps = read_profile_maps(path)
    if not maps:
        sys.exit(f"no profile map in {path}")
    vol = maps[-1]
    phis = range(-(vol.shape[0] // 2), vol.shape[0] // 2 + 1)
    fig, axes = plt.subplots(1, vol.shape[0], figsize=(14, 2.2), sharey=True)
    for ax, sl, phi in zip(axes, vol, phis):
        ax.imshow(sl, cmap="magma", vmin=0, vmax=100, origin="upper")
        ax.set_title(f"φ {phi:+d}")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
