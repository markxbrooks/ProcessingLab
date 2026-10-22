#!/usr/bin/env python3
"""Plot resolution-bin completeness and I/σ from a ProcessingLab log.

Reads a ``RESOLUTION_BINS`` block when the log has one, and the older
XDS-style table (percent signs glued to the numbers) otherwise. The last
table in the file is the one that is plotted.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

COLUMNS = [
    "d_limit_A",
    "observed",
    "unique",
    "possible",
    "completeness_pct",
    "r_obs_pct",
    "r_exp_pct",
    "compared",
    "mean_isigma",
    "rmeas_pct",
    "cc_half_pct",
    "cc_sig",
    "anomal_corr_pct",
    "sigano",
    "nano",
]


def _is_number(token):
    try:
        float(token)
    except ValueError:
        return False
    return True


def _read_marked(lines):
    tables = []
    i = 0
    while i < len(lines):
        if not lines[i].startswith("RESOLUTION_BINS "):
            i += 1
            continue
        cols = lines[i].split("columns=", 1)[1].split(",")
        rows = []
        i += 1
        while i < len(lines) and lines[i] != "RESOLUTION_BINS_END":
            if lines[i].strip() and not lines[i].startswith("#"):
                rows.append(lines[i].split())
            i += 1
        if rows:
            tables.append(pd.DataFrame(rows, columns=cols))
        i += 1
    return tables


def _legacy_frame(raw_rows):
    cleaned = []
    for parts in raw_rows:
        cc = parts[10]
        sig = "1" if cc.endswith("*") else "0"
        bare = [p.rstrip("%*") for p in parts]
        # d observed unique possible compl r_obs r_exp compared I/σ Rmeas CC
        # then Anomal SigAno Nano. cc_sig is not its own column in this form.
        cleaned.append(bare[:11] + [sig] + bare[11:])
    return pd.DataFrame(cleaned, columns=COLUMNS)


def _read_legacy(lines):
    """XDS-style rows: a ``%`` on the completeness field, ``*`` glued to CC."""
    tables = []
    rows = []

    def flush():
        if rows:
            tables.append(_legacy_frame(rows))
            rows.clear()

    for line in lines:
        parts = line.split()
        if (
            len(parts) >= 14
            and (parts[0] == "total" or _is_number(parts[0]))
            and parts[4].endswith("%")
        ):
            rows.append(parts[:14])
            continue
        flush()
    flush()
    return tables


def read_resolution_bins(path):
    """Return every resolution table in ``path``, oldest first.

    :param path: Log file.
    :returns: One dataframe per table. Empty when the log has none.
    """
    lines = Path(path).read_text(errors="replace").splitlines()
    marked = _read_marked(lines)
    if marked:
        return marked
    return _read_legacy(lines)


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


def main():
    path = _resolve_log(sys.argv)
    if not path.is_file():
        sys.exit(f"no log at {path}")
    tables = read_resolution_bins(path)
    if not tables:
        sys.exit(f"no resolution table in {path}")
    df = tables[-1]
    shells = df[df.d_limit_A != "total"].apply(pd.to_numeric)

    fig, ax = plt.subplots()
    ax.plot(shells.d_limit_A, shells.completeness_pct, "o-", color="C0",
            label="completeness %")
    ax.set_xlim(shells.d_limit_A.max(), shells.d_limit_A.min())
    ax.set_xlabel("d (Å)")
    ax.set_ylabel("completeness %")
    ax2 = ax.twinx()
    ax2.plot(shells.d_limit_A, shells.mean_isigma, "s--", color="C1",
             label="I/σ")
    ax2.set_ylabel("I/σ")
    handles = ax.get_lines() + ax2.get_lines()
    ax.legend(handles, [h.get_label() for h in handles], loc="upper right")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
