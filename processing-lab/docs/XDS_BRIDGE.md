# XDS as an optional external bridge

ipmosflm and processing-lab **do not** require [XDS](https://wiki.uni-konstanz.de/xds/index.php/Installation) to build or run. XDS is useful when you need to **convert** an unsupported detector format into miniCBF (or process HDF5 with vendor plugins), then feed the result into MOSFLM’s existing CBFlib path ([CBF_STATUS.md](CBF_STATUS.md)).

## When to use it

| Situation | Suggested workflow |
| --------- | ------------------ |
| Native CBF / miniCBF | Read directly (`DETECTOR CBF` / Pilatus / `cbfwrap_`) |
| Format XDS understands but MOSFLM does not | Convert with XDS `2cbf` (or related tools) → miniCBF → MOSFLM |
| Eiger HDF5 | Prefer MOSFLM `fromhdf5`, **or** convert frames externally (`H5ToXds` / `hdf2mini-cbf`) if needed |
| `*.cbf.gz` | XDS `xds-zcbf` plugin can speed XDS; for MOSFLM, decompress or convert first |

## Install notes (operators only)

Follow the [XDSwiki Installation](https://wiki.uni-konstanz.de/xds/index.php/Installation) page:

- Academic download via Dectris form; commercial users need a subscription.
- macOS Apple Silicon: `XDS-Apple_M1.zip` → `/usr/local/bin` or `$HOME/bin`.
- Quarantine: `sudo xattr -dr com.apple.quarantine XDS-Apple_M1/*` if Gatekeeper blocks binaries.
- HDF5 plugins (Durin / Neggia) and helpers (`generate_XDS.INP`, XDSGUI) are optional extras from the wiki `get_folder_silicon.sh` scripts.
- CCP4 helps XDSSTAT / XDSCONV / POINTLESS from XDSGUI; it is **not** required to compile ipmosflm’s CBF reader.

## Explicit non-goals

- Do **not** vendor XDS binaries into this repository.
- Do **not** add Neggia/Durin as link-time dependencies of `bin/ProcessingLab`,
  `bin/ipmosflm` (oracle), or `processing-lab`.
- Do **not** treat `2cbf` as a runtime requirement for every image open.

## Importing Xia2 `.xinfo` / `XDS.INP`

ProcessingLab can seed CONTROL from Xia2/XDS text files without running XDS:

```bash
./ProcessingLab --xinfo automatic.xinfo --xds-inp path/to/XDS.INP HKLOUT out.mtz script.in
```

| Flag | Source | Maps into CONTROL |
|------|--------|-------------------|
| `--xds-inp` | `XDS.INP` | `WAVELENGTH`, `DISTANCE`, `PIXEL`, `BEAM` (from `ORGX`/`ORGY` x `QX`/`QY` when present), `CELL`, `SYMMETRY`, `DIRECTORY`+`TEMPLATE` (`?` becomes `#`), `PROCESS` (+ `ANGLE` from `OSCILLATION_RANGE`) |
| `--xinfo` | Xia2 `.xinfo` | `WAVELENGTH`, `DIRECTORY`, `IMAGE`/`TEMPLATE`, `PROCESS` from first sweep |

Apply order is **XDS.INP then `.xinfo`**, so local Xia2 paths and `START_END` override stale absolute `NAME_TEMPLATE_OF_DATA_FRAMES` paths.

Orientation is not in `XDS.INP` / `.xinfo`; use `--xparm XPARM.XDS` (or CONTROL `XPARM` / `GXPARM`) to seed `CELL`, `SYMMETRY`, `DISTANCE`, `PIXEL`, `BEAM`, and `MATRIX` when present. `AUTOINDEX` GO writes `XPARM.XDS`; `PROCESS` GO writes `GXPARM.XDS` (disable with `XPARMOUT OFF`). Prefer an index-stage `XDS.INP` (or the image header) for beam centre when the integrate/CORRECT `XDS.INP` omits `ORGX`/`ORGY`.

CORRECT-lite: `POLARIZATION ρ` (default 0.98) in LP; `AIR` path-length absorption; systematic-absence I/σ sanity; `CELLFIX` before PROCESS.

INTEGRATE-lite: `PARTIALSCALE` (I÷FRACTIONCALC), `AUTOSIGMA` (σM/σD from partials), PEAKHALF clamped to half SEPARATION, trusted-mask profile bg.

### IDXREF / AUTOINDEX (COLSPOT coverage)

XDS IDXREF uses multi-image COLSPOT (often thousands of spots). ProcessingLab
`FINDSPOTS` defaults to `MAXSPOTS 1500`. For large cells, **single-frame
`MAXSPOTS ≤ 500` is unsafe** — DPS can lock onto a short-vector / sublattice
(e.g. 76c1b wrongly indexing as ~43×118×60 mC instead of ~62×187×72 mP).

Recommended CONTROL for indexing:

```text
FINDSPOTS THRESHOLD 5 MAXSPOTS 1500 FIND 1 FILE spots.spt MULTIFRAME RANGE 1 50
GO
AUTOINDEX DPS THRESHOLD 15 SEPMIN 6 CLUSTER_RADIUS 2 …
GO
```

`MULTIFRAME` without `RANGE` defaults to a 50-frame wedge from `FIND`.
For CBF/SMV templates, each frame path is resolved separately (HDF5 masters
still use the frame index). Default `MERGE_XY 0` keeps per-frame 3D spots
and **φ-stratifies** `MAXSPOTS` so the brightest frames do not starve later φ.
Set `MERGE_XY` / `MERGE_PHI` only when clustering persistent peaks.
DiffractionLab `build_findspots_commands(..., for_indexing=True)` emits the
same MULTIFRAME pattern with `MAXSPOTS 1500`.

#### DPS + Δ-vector + FFT3D (INDEXFIT)

AUTOINDEX runs **DPS**, a Kabsch-style **difference-vector** search, and a
native **DIALS-style 3D-FFT** lattice search, then ranks candidates with
INDEXFIT (`%indexed` / RMSD HKL). No CCTBX/DIALS link — FFT is pure C++.

| Knob | Default | Meaning |
|------|---------|---------|
| `SEPMIN` | 6 | XDS units of `0.001/WAVE` (Å⁻¹). Auto-capped via `DMAX` so large-cell `\|b*\|` is not excluded. |
| `CLUSTER_RADIUS` | 2 | Same units as SEPMIN. |
| `MININDEXED` | 0.70 | XDS-style soft gate (≥70% indexed). Below → `INDEXFIT: INSUFFICIENT` warning. |
| `DELTA` / `NODELTA` | on | Enable / disable the Δ-vector candidate path. |
| `FFT` / `FFT3D` / `NOFFT` | on | Enable / disable the native 3D-FFT path. |
| `DPSONLY` | — | DPS only (disables Δ-vector and FFT). |
| `FFTONLY` | — | FFT only (debug; skips DPS and Δ-vector). |

Look for `DIFFVECTOR:` / `FFT3D:` / `INDEXFIT: method=dps|delta|fft3d` in the
log. Explicit `SYSTEM` / `CHARACTER` / `SOLUTION` (e.g. insulin `SOLUTION 12`
→ cI) keeps the DPS Bravais path and skips Δ-vector / FFT override.

FFT needs denser φ coverage than classic single-frame DPS. The native FFT
keeps a **fixed real-space window** up to ``DMAX`` (only low-|q| spots enter
the histogram) so cell edges are not crushed to ~10 Å junk. Basis search
follows DIALS: seed the **three shortest non-collinear** peaks, inject
``(u±v)/2`` candidates, then apply a ModelRankFilter-style doubled-cell /
volume (≤1.25×) / ``n_indexed`` (≥0.9×) gate so face-diagonal **√2
supercells** (~84 Å instead of ~58 Å for thaumatin) lose to the primitive.
For large cells:

```text
FINDSPOTS THRESHOLD 5 MAXSPOTS 20000 FIND 1 FILE spots.spt MULTIFRAME RANGE 1 900
GO
AUTOINDEX DPS THRESHOLD 12 DMAX 200 FFT …
GO
```

BESTSOLAR-lite among LASTGOOD still prefers `mP` over a poorly indexing `mC`.
If AUTOINDEX still fails on short φ wedges, seed from Xia2/XDS with `--xparm`
pointing at an **index-stage** `XPARM.XDS` (not a prior bad ProcessingLab
`GXPARM.XDS`). Harden CI with `MOSFLM_76C1B_REQUIRE_CELL=1` once spot-only
recovery is solid.

Optional multi-image COLSPOT: `FINDSPOTS MULTIFRAME RANGE first last …`. Multi-frame PROCESS runs periodic post-refinement (`POSTREFINE PERIOD n`, default every 25 frames).

### NPROC (OpenMP frame parallelism)

When built with OpenMP (`brew install libomp` on macOS), CONTROL `NPROC n`
(aliases `THREADS` / `NTHREADS`) parallelizes:

| Stage | How |
|-------|-----|
| `FINDSPOTS MULTIFRAME` | `#pragma omp parallel for` over frames; quota / merge stays serial |
| `PROCESS first TO last` | Parallel inside each `POSTREFINE PERIOD` chunk; postrefine + `AUTOSIGMA` stay serial barriers |

`NPROC OFF` / `1` forces serial. Default is `PROCESSINGLAB_NPROC` env, else
hardware max (clamped to 32). Serial builds ignore `NPROC` (log shows
`(serial build)`). Example:

```text
NPROC 8
FINDSPOTS … MULTIFRAME RANGE 1 50
…
PROCESS 1 TO 900
GO
```
