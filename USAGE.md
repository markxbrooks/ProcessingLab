# ProcessingLab usage

`ProcessingLab` is a MOSFLM-style CONTROL driver. It reads images and a
keyword script (or Xia2/XDS text files that are turned into those keywords),
then writes unmerged MTZ plus XDS-style orientation snapshots.

It does **not** run XDS. `--xds-inp` and `--xinfo` only parse those files and
seed CONTROL. Native XDS binaries are optional and unused at runtime.

Build the binary with `make cpp`. The install path is `bin/ProcessingLab`.

```bash
./bin/ProcessingLab --help
```

## Quick start: `ProcessingLab run`

```bash
ProcessingLab run --image /data/sweep_master.h5                 # whole sweep
ProcessingLab run --image /data/sweep_0001.cbf --range 1:180 \
    --xds-inp /data/XDS.INP --output-dir results/lysozyme
```

`run` finds spots, autoindexes and integrates each sweep without a CONTROL
script. The image headers supply the experiment metadata; XDS.INP and
`.xinfo` override or add to it.

| Option | Role |
| ------ | ---- |
| `--image PATH` | HDF5/NeXus master, or one frame of a CBF/SMV sweep (its numbered siblings are the sweep). Repeatable: each image is one sweep, processed in turn into `sweep_1/`, `sweep_2/`, … Optional if `--xinfo` / `--xds-inp` names the images. |
| `--range FIRST:LAST` | Inclusive frame range. Default: every frame found. |
| `--output-dir DIR` | Directory for all outputs. It must not exist yet; runs never reuse one. Default: a new `processing_YYYYMMDD_HHMMSS` (with `_2`, `_3`, … if taken) in the current directory. |
| `--xds-inp PATH` | Import XDS.INP. |
| `--xinfo PATH` | Import Xia2 `.xinfo`. |
| `--sweep NAME` | Which `.xinfo` sweep, when it has several (required then). |
| `--help`, `--version` | Usage, build version. |

**Metadata precedence**, per field, lowest to highest: image header, XDS.INP,
`.xinfo`, command line. A source only overrides the fields it sets (an
`.xinfo` without detector geometry keeps the header's). `process.log` starts
with the resolved values and their sources, and warns when sources disagree
on a value (the higher-priority one is used):

```text
 Experiment metadata
  Image source:  /data/sweep_0001.cbf  [CLI]
  Wavelength:    0.97949 Å  [XDS.INP]
  Distance:      384.99 mm  [XDS.INP]
  Beam:          214.28 207.96 mm  [image header]
  Frame range:   1–180  [CLI]  (900 on disk)
  *** WARNING: Distance: XDS.INP 384.99 mm overrides image header 395.00 mm
```

Every output goes in the output directory: `input.json`, `process.json`,
`process.md`, `process.log`, `process.mtz`,
`spots.spt`, `XPARM.XDS`, `GXPARM.XDS`, `ProcessingLabIndex.txt`. Spots
for indexing come from the first 50 frames of the range. With a known cell
(from XDS.INP), `DMAX` is set to 1.7× its longest edge so a long axis is not
lost in autoindexing. For anything else (fixed matrix, resolution limits,
`NPROC`, …) use a CONTROL script as below.

## Invocation

```text
ProcessingLab run [options]          (see above)
ProcessingLab [HKLOUT out.mtz] [--image img]
              [--xinfo file.xinfo] [--xds-inp XDS.INP]
              [--xparm XPARM.XDS] [--output-dir DIR] [script.in]
```

CONTROL is read from `script.in` if given, otherwise from stdin. `GO` / `RUN`
runs any pending `FINDSPOTS`, `AUTOINDEX`, or `PROCESS`. `EXIT` / `QUIT` ends
the session.

| Argument | Role |
| -------- | ---- |
| `HKLOUT path` | Unmerged MTZ path. Same as the CONTROL keyword. Default `process.mtz`. |
| `--image path` | Single image or HDF5/NeXus master. Overrides `IMAGE` unless `TEMPLATE` is set for multi-frame CBF/SMV. |
| `--xds-inp path` | Parse `XDS.INP` into CONTROL (see below). |
| `--xinfo path` | Parse Xia2 `.xinfo` into CONTROL (first sweep and the wavelength it names). `run --sweep` selects another. |
| `--xparm path` | Parse `XPARM.XDS` / `GXPARM.XDS` into `CELL`, `SYMMETRY`, geometry, and `MATRIX`. Alias `--gxparm`. |
| `--output-dir DIR` | New directory (must not exist) for `process.log`, the MTZ, spot lists, `XPARM.XDS`, `GXPARM.XDS` and the run records (`input.json`, `process.json`, `process.md`). Relative paths in the script (`DIRECTORY`, `IMAGE`) still resolve against the current directory. |
| `script.in` | CONTROL file. Lines starting with `!` are comments. |

Apply order before the script: **`--xds-inp`, then `--xinfo`, then `--xparm`**.
Later sources overwrite earlier ones for the same field. Script keywords
overwrite the imports.

If no `--image`, `IMAGE`, `TEMPLATE`, `--xinfo`, or `--xds-inp` supplies a
path, ProcessingLab falls back to `$MOSFLM_IMAGE` or
`$CCP4/examples/data/insulin_1_001.img` when that file exists.

OpenMP workers: CONTROL `NPROC n` (aliases `THREADS` / `NTHREADS`), or the
`PROCESSINGLAB_NPROC` environment variable. `NPROC OFF` / `1` is serial.

## Images

| Format | How to point at it |
| ------ | ------------------ |
| HDF5 / NeXus Eiger master (`*_master.h5`, `.nxs`) | `--image master.h5`. Frame numbers are indices in the master. Keep `data_*.h5` beside it. |
| miniCBF / CBF | `DIRECTORY` + `TEMPLATE` with `#` placeholders, e.g. `sweep_####.cbf`. |
| SMV (ADSC) | Same `DIRECTORY` / `TEMPLATE` / `--image` pattern. Read as MOSFLM does: `XDET` runs along the file's slow axis (DIALS `XDET` is the fast axis). |

For multi-frame CBF/SMV (`FINDSPOTS … MULTIFRAME`), each frame is read from
the `--image` file's sibling with that frame number (same zero padding), else
from `TEMPLATE`; it is an error if neither exists. HDF5 masters are fine as a
single `--image`.

Geometry not set in CONTROL is filled from the image header on first read
(distance, wavelength, beam, pixel size).

## CONTROL script

A typical index-and-integrate run:

```text
NPROC 8
DIRECTORY /data/sweep
TEMPLATE sweep_####.cbf
WAVELENGTH 0.9795
DISTANCE 200.0
PIXEL 0.075 0.075
BEAM 166.16 172.97
FINDSPOTS THRESHOLD 5 MAXSPOTS 20000 FIND 2 FILE spots.spt MULTIFRAME RANGE 2 180
GO
AUTOINDEX DPS THRESHOLD 12 HKLDEV 0.3 NVECT 40 NVREF 40 SEPMIN 6 CLUSTER_RADIUS 2 FILE spots.spt
GO
PROCESS 2 TO 751
GO
EXIT
```

`GO` after `FINDSPOTS` writes the spot file; `GO` after `AUTOINDEX` indexes
it and writes `XPARM.XDS`; `GO` after `PROCESS` integrates and writes MTZ /
`GXPARM.XDS`.

### Geometry and crystal

| Keyword | Meaning |
| ------- | ------- |
| `DISTANCE` | Crystal-to-detector (mm). |
| `WAVELENGTH` / `WAVE` | Å. |
| `PIXEL px [py]` | Pixel size (mm). |
| `BEAM x y` | Beam centre (mm). |
| `CELL a b c α β γ` | Direct cell (Å, degrees). |
| `SYMMETRY n` | Space-group number (1–230). |
| `MATRIX` then 9 numbers | MOSFLM orientation matrix [A] (lab X,Y,Z rows). |
| `CELLFIX ON\|OFF\|HARD\|SOFT [s]` | Project the cell onto the space-group metric before PROCESS (default ON, hard). |
| `MOSAIC` | Mosaic spread (degrees). Default 0.15. |
| `DIVERGENCE` / `DIVH` `DIVV` | Beam divergence (degrees). |
| `YSCALE` `TILT` `TWIST` `ROTMAT2` | Detector fiducials (insulin-style). |

### Images and scan

| Keyword | Meaning |
| ------- | ------- |
| `DIRECTORY` | Working directory for relative `IMAGE` / `TEMPLATE` / spot files. |
| `IMAGE` | Single-frame path. |
| `TEMPLATE` | Multi-frame basename with `#` digits (`76c1b_1_####.cbf`). |
| `PROCESS first TO last [START φ] [ANGLE Δφ]` | Frames to integrate. `START`/`ANGLE` override header φ. |
| `DMIN` / `DMAX` / `RESOLUTION` | Predictor resolution limits (Å). Unset `DMIN` uses the farthest inset detector corner. An explicit value overrides that. |
| `REFINE ALL` or `REFINE ORIENTATION CELL DISTANCE BEAM TILT` (any subset) | Parameters positional refinement may change during PROCESS. Default `ALL`; `--xparm` sets `ORIENTATION` (as XDS INTEGRATE; the imported cell and detector are kept). |
| `INVERTX ON\|OFF` | Handedness of a supplied `MATRIX`: `ON` (default) for MOSFLM/ADSC matrices, `OFF` for ProcessingLab or XDS-frame matrices (`--xparm` sets it). |
| `PHI_OFFSET` | Added to every frame’s φ. |

### Spot search and indexing

| Keyword | Meaning |
| ------- | ------- |
| `FINDSPOTS THRESHOLD σ MAXSPOTS n FIND f FILE path.spt [MULTIFRAME RANGE a b]` | COLSPOT analogue. Default `MAXSPOTS` 1500; use much more for large cells. `MULTIFRAME` without `RANGE` is a 50-frame wedge from `FIND`. |
| `AUTOINDEX DPS … FILE path.spt` | Index the `.spt`. If `FILE` is omitted, the previous `FINDSPOTS` path is used. |

Useful `AUTOINDEX` knobs: `THRESHOLD`, `HKLDEV`, `NVECT`, `NVREF`, `SEPMIN`,
`CLUSTER_RADIUS`, `MININDEXED` (default 0.70), `SYSTEM` / `CHARACTER` /
`SOLUTION`, `DELTA`/`NODELTA`, `FFT`/`NOFFT`, `DPSONLY`, `IDXREF`/`NOIDXREF`.
Indexing details live in [`processing-lab/docs/XDS_BRIDGE.md`](processing-lab/docs/XDS_BRIDGE.md).

### Integration

| Keyword | Meaning |
| ------- | ------- |
| `SEPARATION x [y]` | Measurement-box size (mm). |
| `PEAKHALF n` | Peak half-width (pixels); clamped to half `SEPARATION`. |
| `MINFRAC` | Drop reflections with `FRACTIONCALC` below this (default 0.5). |
| `PARTIALSCALE ON\|OFF` | Scale `I`/`σ` by `1/FRACTIONCALC` (default ON). |
| `POLARIZATION ρ` | LP polarisation fraction (default 0.98). |
| `AIR c` | Air absorption (mm⁻¹). |
| `AUTOSIGMA ON\|OFF` | Estimate mosaic / divergence from partials (default ON). |
| `POSTREFINE [PERIOD n] \| OFF` | Misset refine every *n* frames (default period 25). |
| `WIDEN_EWALD ON\|OFF` | Widen Ewald with mosaic (PROCESS default OFF). |
| `NMASKS nx ny` | Standard-profile detector bins. |

`XPARM path` / `GXPARM path` in the script is the same import as `--xparm`.
`XPARMOUT path` writes a snapshot immediately; `XPARMOUT OFF` disables the
automatic `GXPARM.XDS` after PROCESS.

## `XDS.INP`

`--xds-inp` reads an XDS-style `KEY=value` file (`!` comments). Recognised
keys:

| `XDS.INP` | CONTROL |
| --------- | ------- |
| `X-RAY_WAVELENGTH` | `WAVELENGTH` |
| `DETECTOR_DISTANCE` | `DISTANCE` |
| `QX` / `QY` | `PIXEL` |
| `ORGX` / `ORGY` plus `QX`/`QY` | `BEAM` in mm (`ORG · pixel size`) |
| `UNIT_CELL_CONSTANTS` | `CELL` |
| `SPACE_GROUP_NUMBER` | `SYMMETRY` |
| `DATA_RANGE` | `PROCESS first TO last` |
| `OSCILLATION_RANGE` | `PROCESS … ANGLE` |
| `NAME_TEMPLATE_OF_DATA_FRAMES` | `DIRECTORY` + `TEMPLATE` (`?` → `#`) |

`ORGX`/`ORGY` without pixel size do **not** set `BEAM`. Integrate/CORRECT
`XDS.INP` files often omit them; prefer an IDXREF-stage `XDS.INP` or the
image header. Orientation is never taken from `XDS.INP` — use `--xparm`.

Example (from [`processing-lab/testdata/import/sample_xds.inp`](processing-lab/testdata/import/sample_xds.inp)):

```text
ORGX=1245.82 ORGY=1209.07
SPACE_GROUP_NUMBER=3
UNIT_CELL_CONSTANTS= 61.70 185.40 71.90 90.00 100.20 90.00
NX=2463 NY=2527 QX=0.172 QY=0.172
DETECTOR_DISTANCE=385.190
OSCILLATION_RANGE=0.20
X-RAY_WAVELENGTH=0.979490
NAME_TEMPLATE_OF_DATA_FRAMES=/data/results/old/sweep_1_0???.cbf
DATA_RANGE=1 900
```

becomes `PIXEL 0.172 0.172`, `BEAM 214.28 207.96`, `TEMPLATE sweep_1_0###.cbf`,
`PROCESS 1 TO 900 ANGLE 0.2`, and so on.

You still need `GO` (in a script or on stdin) to run PROCESS. A minimal
follow-up when geometry is already imported:

```bash
./bin/ProcessingLab HKLOUT out.mtz --xds-inp XDS.INP <<'EOF'
GO
EXIT
EOF
```

To index first, ignore imported `CELL`/`SYMMETRY`/`PROCESS` by running
`FINDSPOTS` / `AUTOINDEX` in the script (those keywords replace the pending
`PROCESS` only when you issue a new `PROCESS` after indexing — AUTOINDEX
fills `CELL`/`MATRIX`/`SYMMETRY` itself). Typical pattern:

```bash
./bin/ProcessingLab HKLOUT out.mtz \
  --xds-inp XDS.INP --xinfo automatic.xinfo <<'EOF'
FINDSPOTS THRESHOLD 5 MAXSPOTS 1500 FIND 1 FILE spots.spt MULTIFRAME RANGE 1 50
GO
AUTOINDEX DPS THRESHOLD 15 SEPMIN 6 CLUSTER_RADIUS 2
GO
PROCESS 1 TO 900
GO
EXIT
EOF
```

## Xia2 `.xinfo`

`--xinfo` reads the **first** `BEGIN WAVELENGTH` block and the **first**
`BEGIN SWEEP` block. `!` starts a comment; `#` in paths is kept.

| `.xinfo` | CONTROL |
| -------- | ------- |
| `WAVELENGTH` inside the wavelength block | `WAVELENGTH` |
| Sweep `DIRECTORY` | `DIRECTORY` |
| Sweep `IMAGE` | `IMAGE`, or `TEMPLATE` when `START_END` spans more than one frame (digits before the extension become `#`) |
| `START_END first last` | `PROCESS first TO last` |

Example (from [`processing-lab/testdata/import/sample.xinfo`](processing-lab/testdata/import/sample.xinfo)):

```text
BEGIN WAVELENGTH NATIVE
WAVELENGTH 0.979490
END WAVELENGTH NATIVE

BEGIN SWEEP SWEEP1
WAVELENGTH NATIVE
DIRECTORY /data/sweep
IMAGE sweep_1_0001.cbf
START_END 1 900
END SWEEP SWEEP1
```

Applied **after** `--xds-inp`, so Xia2 `DIRECTORY` / `START_END` / `IMAGE`
override stale absolute `NAME_TEMPLATE_OF_DATA_FRAMES` paths from XDS.

`.xinfo` has no detector distance, pixel size, beam, cell, or matrix. Combine
with `--xds-inp` and/or `--xparm`, or let the image header fill geometry.

## `XPARM.XDS` / `GXPARM.XDS`

`--xparm` (or CONTROL `XPARM` / `GXPARM`) seeds:

- `CELL`, `SYMMETRY`, `WAVELENGTH`, `DISTANCE`, `PIXEL`, `BEAM` (from `ORGX`/`ORGY` × pixel size)
- `MATRIX` when the file contains an orientation, with `INVERTX OFF` and
  `REFINE ORIENTATION`

Use an **index-stage** `XPARM.XDS` (Xia2 `…/index/XPARM.XDS` or ProcessingLab’s
own post-AUTOINDEX file). A previous bad `GXPARM.XDS` will lock PROCESS onto
the wrong lattice.

## Outputs

With `--output-dir` (always, for `run`) every relative path is in the output
directory. Otherwise the log and MTZ are relative to the current directory,
and spot lists, `XPARM.XDS` and `GXPARM.XDS` to `DIRECTORY` (default `.`).

| File | When | Contents |
| ---- | ---- | -------- |
| `input.json` | start (with an output directory) | Resolved inputs: images, imports, each experiment value as `{value, unit, source, overrides}`, effective FINDSPOTS / AUTOINDEX / integration / refinement settings, the CONTROL run. |
| `process.json` | created at start, updated at every `GO` | `status` (`running` / `completed` / `failed`), `stages` with status and timing, `indexing.selected_solution`, `refinement` (after PROCESS), `integration`, `quality` (unscaled shells + total), `warnings`, `error`, `outputs`. Written atomically. |
| `process.md` | with `process.json` | Report rendered from the two records. |
| `process.log` | always (unless `--no-log`) | Everything printed to the terminal. |
| `HKLOUT` (`process.mtz` if unset) | PROCESS `GO` | Unmerged MOSFLM-style MTZ. |
| `spots.spt` (or `FINDSPOTS FILE`) | FINDSPOTS `GO` | MOSFLM spot list: header (nfast, nslow, pixel, yscale, omega; invflag; beam mm) then `x_mm y_mm φ_start φ_end I σ`. |
| `XPARM.XDS` | successful AUTOINDEX | Bridge-format orientation + cell + detector (written next to `DIRECTORY`). |
| `GXPARM.XDS` | PROCESS `GO` (unless `XPARMOUT OFF`) | Same layout after integration. Override the path with `XPARMOUT`. |
| `ProcessingLabIndex.txt` | successful AUTOINDEX | Sidecar `CELL` / `MATRIX` / `SYSTEM` / `CHARACTER` / `INDEXED` / `RMSD_HKL`. |
| stdout | always | Banner, import lines, COLSPOT/INDEXFIT/IDXREF/PROCESS log, `PROCESSINGLAB_*` dump. |

### MTZ columns

PROCESS writes Gemmi MTZ titled `processing-lab PROCESS` with space group and
cell from CONTROL (after `CELLFIX`). Columns match MOSFLM `STARTMTZ`:

`H K L` `M/ISYM` `BATCH` `I` `SIGI` `IPR` `SIGIPR` `FRACTIONCALC`
`XDET` `YDET` `ROT` `WIDTH` `LP` `MPART` `FLAG` `BGPKRATIOS`

`BATCH` is the image number. Batch headers are present (cell + wavelength)
so POINTLESS/AIMLESS treat the file as unmerged. `I` is the working
intensity (profile-fitted when the standard-profile fitter succeeds,
otherwise summation); `IPR` is the profile estimate. Reflections with
`FRACTIONCALC < MINFRAC` are omitted.

### Machine-readable AUTOINDEX lines

```text
PROCESSINGLAB_CELL a b c α β γ
PROCESSINGLAB_MATRIX a11 … a33
PROCESSINGLAB_SYSTEM tP
PROCESSINGLAB_CHARACTER 10
PROCESSINGLAB_NSPOTS 5000
PROCESSINGLAB_INDEXED 96
PROCESSINGLAB_RMSD_HKL 0.034
```

## Worked examples

HDF5 master, spots from scratch:

```bash
./bin/ProcessingLab HKLOUT process.mtz --image ~/data/thaumatin_9_1_master.h5 <<'EOF'
FINDSPOTS THRESHOLD 5 MAXSPOTS 20000 FIND 2 FILE spots.spt MULTIFRAME RANGE 2 180
GO
AUTOINDEX DPS THRESHOLD 12 FILE spots.spt
GO
PROCESS 2 TO 751
GO
EXIT
EOF
```

Xia2 + XDS files, then integrate the imported range:

```bash
./bin/ProcessingLab HKLOUT out.mtz \
  --xds-inp /path/to/XDS.INP \
  --xinfo /path/to/automatic.xinfo \
  --xparm /path/to/index/XPARM.XDS <<'EOF'
GO
EXIT
EOF
```

Known matrix, no indexing (insulin-style):

```bash
./bin/ProcessingLab HKLOUT /tmp/out.mtz --image insulin_1_001.img \
  processing-lab/testdata/images/known_indexed/cpp_process.in
```
