# ProcessingLab usage

`ProcessingLab` is a MOSFLM-style CONTROL driver. It reads images and a
keyword script (or Xia2/XDS text files that are turned into those keywords),
then writes unmerged MTZ plus XDS-style orientation snapshots.

It does **not** run XDS. `--xds-inp` and `--xinfo` only parse those files and
seed CONTROL. Native XDS binaries are optional and unused at runtime.

## Build and install

```bash
./configure --prefix=$HOME/.local   # default prefix: the CCP4 tree ($CCP4)
make                                # → bin/ProcessingLab in the source tree
make install                        # → $(prefix)/bin/ProcessingLab
./bin/ProcessingLab --help
```

`configure` needs `$CCP4` set, or `--with-ccp4-prefix=/path/to/ccp4`.
`make` skips DPS autoindexing, with a warning, when gfortran is missing.
`make install` only copies the binary `make` built. Use `sudo` for
`make install` alone, and only when the prefix needs it. `make uninstall`
removes it again.

## Quick start: `ProcessingLab run`

```bash
ProcessingLab run --image /data/sweep_master.h5                 # whole sweep
ProcessingLab run --image /data/sweep_0001.cbf --range 1:180 \
    --xds-inp /data/XDS.INP --output-dir results/lysozyme
```

`run` finds spots, autoindexes, integrates and runs `CORRECT` on each sweep
without a CONTROL script. The image headers supply the experiment metadata; XDS.INP and
`.xinfo` override or add to it.

| Option | Role |
| ------ | ---- |
| `--image PATH` | HDF5/NeXus master, or one frame of a CBF/SMV sweep (its numbered siblings are the sweep). Repeatable: each image is one sweep, processed in turn into `sweep_1/`, `sweep_2/`, … Optional if `--xinfo` / `--xds-inp` names the images. |
| `--range FIRST:LAST` | Inclusive frame range. Default: every frame found. |
| `--output-dir DIR` | Directory for all outputs. It must not exist yet; runs never reuse one. Default: a new `processing_YYYYMMDD_HHMMSS` (with `_2`, `_3`, … if taken) in the current directory. |
| `--xds-inp PATH` | Import XDS.INP. |
| `--xinfo PATH` | Import Xia2 `.xinfo`. |
| `--sweep NAME` | Which `.xinfo` sweep, when it has several (required then). |
| `--xparm PATH` | Import `XPARM.XDS` / `GXPARM.XDS`: cell, space group, geometry and orientation. With an orientation, autoindexing is skipped and only the orientation is refined. |
| `--from RUN` | Re-process from a completed run: its directory or `process.json`. See [Re-processing](#re-processing-with-refined-parameters-run---from). |
| `--spacegroup N` | Space group number or symbol (`92`, `P41212`), over every import and `--from`. See [Overriding](#overriding-the-space-group-and-cell-run---spacegroup-run---cell). |
| `--cell A B C α β γ` | Unit cell, over every import and `--from`; must fit the space group. |
| `--refine LIST` | What PROCESS refines: `orientation`, `cell`, `distance`, `beam`, `tilt_twist`, `mosaicity` (re-estimated), or `all` / `none`; the rest are held fixed. Default: `orientation,mosaicity` with `--from` / `--xparm`, otherwise `all`. |
| `--untrusted SHAPE V…` | Detector region to ignore, in XDS pixels (repeatable): `RECTANGLE x1 x2 y1 y2`, `ELLIPSE x1 x2 y1 y2` (its bounding box) or `QUADRILATERAL x1 y1 x2 y2 x3 y3 x4 y4`. Added to any `UNTRUSTED_*=` of `--xds-inp`, and written to `xds/XDS.INP`. See [Untrusted pixels](#untrusted-pixels-defpix). |
| `--exclude-resolution D1 D2` | Resolution range (Å) to ignore (repeatable). |
| `--exclude-ice` | Ignore the usual ice-ring ranges. |
| `--run-mode MODE` | Where the parameters come from, recorded in the log, `input.json` and `process.md`: `processinglab` (`RunMode.PROCESSINGLAB`, the default: image headers, overridden by any imports), `xds` (`RunMode.XDS`: needs `--xds-inp`, no `--xinfo`), `xia2` (`RunMode.XIA2`: needs `--xinfo`, no `--xds-inp`). `mosflm` (`RunMode.MOSFLM`) is a CONTROL script run (see below). |
| `--help`, `--version` | Usage, build version. |

**Metadata precedence**, per field, lowest to highest: image header, XDS.INP,
`.xinfo`, XPARM (`--xparm`, or the refined values of `--from`), command line. A source only overrides the fields it sets (an
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
`spots.spt`, `XPARM.XDS`, `GXPARM.XDS`, `index.log`, and the
[repeat inputs](#repeat-inputs-xds-xia2-mosflm) `xds/`, `xia2/` and
`mosflm/`. Spots
for indexing come from the first 50 frames of the range. With a known cell
(from XDS.INP), `DMAX` is set to 1.7× its longest edge so a long axis is not
lost in autoindexing. For anything else (fixed matrix, resolution limits,
`NPROC`, …) use a CONTROL script as below.

## Repeat inputs (`xds/`, `xia2/`, `mosflm/`)

Every run, and every CONTROL script run with `--output-dir`, starts by
writing its resolved experiment as input for the other run modes:

```text
RUN/xds/XDS.INP            run --run-mode xds --xds-inp RUN/xds/XDS.INP
RUN/xds/XPARM.XDS          copied after AUTOINDEX (for --xparm)
RUN/xds/GXPARM.XDS         copied after PROCESS (for --xparm)
RUN/xia2/automatic.xinfo   run --run-mode xia2 --xinfo RUN/xia2/automatic.xinfo
RUN/xia2/run_xia2.sh       xia2 xinfo=automatic.xinfo
RUN/mosflm/mosflm.in       --run-mode mosflm --output-dir NEW RUN/mosflm/mosflm.in
```

| File | Contents |
| ---- | -------- |
| `XDS.INP` | `JOB`, `DETECTOR` / `OVERLOAD` (from the header id and saturation value), `NAME_TEMPLATE_OF_DATA_FRAMES` (`?` for frame digits; an Eiger master as `x_??????.h5`), `DATA_RANGE`, `SPOT_RANGE` (first 50 frames, as `run`), `BACKGROUND_RANGE`, `OSCILLATION_RANGE`, `STARTING_ANGLE`, wavelength, distance, `NX NY QX QY`, `ORGX ORGY`, `ROTATION_AXIS` (from the header axis), beam and detector directions, polarisation, and `SPACE_GROUP_NUMBER` / `UNIT_CELL_CONSTANTS` when known. |
| `automatic.xinfo` | One crystal and sweep: `USER_SPACEGROUP` / `USER_CELL` when known, the wavelength, `DIRECTORY`, `IMAGE`, `START_END`, and `BEAM` / `DISTANCE` only where the run overrode the image header. |
| `run_xia2.sh` | `xia2 xinfo=automatic.xinfo`, run from `xia2/`; extra arguments are passed on. |
| `mosflm.in` | The run's CONTROL, self-contained (`IMAGE` names an HDF5 master) and without `HKLOUT`, so outputs go to the new `--output-dir`. |

Read back by ProcessingLab, each reproduces the run it came from. The
`.xinfo` reader understands xia2's `USER_CELL`, `USER_SPACEGROUP`, `BEAM`
and `DISTANCE`, and the XDS.INP reader maps `x_??????.h5` to `x_master.h5`.
The files hold the values at the start of the run. For the refined values
use `run --from`; the copied `XPARM.XDS` / `GXPARM.XDS` are in
ProcessingLab's own layout, for `--xparm`.

## Re-processing with refined parameters (`run --from`)

```bash
ProcessingLab run --from processing_20261022_235027 --output-dir rerun
ProcessingLab run --from processing_20261022_235027/process.json --range 1:450
```

`--from` takes a completed run directory, or its `process.json`, and
processes again from that run's final model:

| From | Values |
| ---- | ------ |
| `process.json` `experiment` | Cell, space group, wavelength, distance, beam, pixel size, tilt/twist, mosaicity, orientation matrix |
| `input.json` | First image and frame range, unless `--image` / `--range` are given |

The model is checked before anything runs. A matrix that implies a cell other
than the run's own is refused, and so is a cell its space group forbids, a
wavelength or pixel size that disagrees with the run's `input.json`, or
missing geometry. Runs from before schema 2 (no `experiment` block) use
their `GXPARM.XDS` and `refinement.tilt_twist_mrad` instead.

Each value is tagged `[previous run]` in `process.log` (`previous_run` in
`input.json`, which also records the run under `inputs.imports`). With the
orientation known, `AUTOINDEX` is skipped, and there is no `XPARM.XDS` or
`index.log` in the output. `FINDSPOTS` still runs, to size the measurement
box. The refinement policy is printed with the metadata:

```text
  Mosaicity:     0.101°  [previous run]  (re-estimated)
  Refine:        orientation, mosaicity; fixed: cell, distance, beam, tilt_twist
```

`--refine` changes it. `--refine orientation,cell` also refines the cell,
and `--refine orientation` holds the mosaicity at its previous value. The
policy is recorded in `input.json` (`experiment.refinement_policy`) and what
was done in `process.json` (`refinement.parameters`). The full contract is
in [processing-lab/docs/RUN_RECORDS.md](processing-lab/docs/RUN_RECORDS.md#second-round-processing-run---from).

```text
  Cell:          57.94 57.94 150.24 90.00 90.00 90.00  [previous run]
  Tilt/twist:    9.39 -8.43 mrad  [previous run]
  Orientation:   refined matrix, no autoindexing  [previous run]
```

Refused, before anything is created: a run whose `process.json` is
`failed` or `running` or has no `refinement`, a multi-sweep run directory
(give its `sweep_N`), and `--from` together with `--xparm`.
`--run-mode xds` / `xia2` refuse both options, as they name a single
parameter source.

`run --xparm PATH` imports an `XPARM.XDS` / `GXPARM.XDS` the same way, at
the same precedence, without the tilt/twist from `process.json`.

## Overriding the space group and cell (`run --spacegroup`, `run --cell`)

```bash
ProcessingLab run --from RUN --spacegroup 92
ProcessingLab run --from RUN --spacegroup P41212 --cell 57.898 57.898 150.952 90 90 90
ProcessingLab run --image /data/x_0001.cbf --spacegroup 19 --cell "62.0,186.9,72.4,90,90,90"
```

| Option | Values | Effect |
| ------ | ------ | ------ |
| `--spacegroup` (alias `--symmetry`) | Number or symbol | `SYMMETRY`, tagged `[CLI]`, over header, imports and `--from`. |
| `--cell` | Six numbers (Å, °), separate or one comma/space list | `CELL`, tagged `[CLI]`. With an orientation from `--from` / `--xparm`, the matrix is rebuilt as U·B(cell): the same crystal orientation (handedness kept) with the new cell. |

The cell must fit the space group's crystal system (monoclinic α = γ = 90°;
orthorhombic all 90°; tetragonal also a = b; hexagonal/trigonal a = b,
γ = 120°; rhombohedral R a = b = c, α = β = γ; cubic a = b = c). Otherwise
the run is refused with the rule it broke. With `--from` only the
orientation is refined, so a given cell stays fixed; if it doesn't match the
data at the refined distance, the positional residual in `process.log`
rises. Changing only the space group within one lattice doesn't alter the
predictions.

## Run modes (`--run-mode`)

Every run has a run mode: where its parameters come from. It is written to
`process.log`, to `input.json` (`run_mode`, `parameter_source`) and to
`process.md`. Inputs that contradict the mode are refused before any
processing.

| `--run-mode` | Shown as | Parameters from | Requires |
| ------------ | -------- | --------------- | -------- |
| `processinglab` | `RunMode.PROCESSINGLAB` | Image headers, overridden per field by any `--xds-inp` / `--xinfo` / `--xparm` / `--from` | Nothing. Default for `run`. |
| `xds` | `RunMode.XDS` | `XDS.INP`; headers fill unset fields | `run` with `--xds-inp`; no `--xinfo` / `--sweep` / `--xparm` / `--from` |
| `xia2` | `RunMode.XIA2` | Xia2 `.xinfo`; headers fill unset fields | `run` with `--xinfo`; no `--xds-inp` / `--xparm` / `--from` |
| `mosflm` | `RunMode.MOSFLM` | MOSFLM CONTROL keywords | A CONTROL script (no `run`). Default for scripts. |

Values are case-insensitive.

```bash
ProcessingLab run --run-mode xds --xds-inp /data/XDS.INP --range 1:180
ProcessingLab run --run-mode xia2 --xinfo automatic.xinfo --sweep SWEEP1
ProcessingLab --run-mode mosflm --image /data/sweep_0001.cbf script.in
```

```text
 Run mode:        RunMode.XDS  (parameters from XDS.INP (--xds-inp); image headers fill unset fields)
```

Mismatches fail with what to change, for example:

```text
ProcessingLab run: RunMode.XDS needs --xds-inp XDS.INP
ProcessingLab run: RunMode.XIA2 takes parameters from the .xinfo only; drop --xds-inp or use --run-mode xds
ProcessingLab: RunMode.XDS is for 'ProcessingLab run'; a CONTROL script is RunMode.MOSFLM
```

## Invocation

```text
ProcessingLab run [options]          (see above)
ProcessingLab [HKLOUT out.mtz] [--image img] [--run-mode mosflm]
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
| `--run-mode mosflm` | Optional: a script run is `RunMode.MOSFLM` (parameters from MOSFLM CONTROL keywords), recorded in the log and `input.json`. Other modes are refused for scripts. |
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
| HDF5 / NeXus Eiger master (`*_master.h5`, `.nxs`) | `--image master.h5`. Frame numbers are indices in the master. Keep `data_*.h5` beside it. Reads DECTRIS `data_NNNNNN` blocks or a single NXmx `data` dataset (virtual datasets too), applies `pixel_mask`, and scales pixel size and distance by their `units`. Reader adapted in part from Diamond's durin ([third_party/durin/](third_party/durin/README.md)). |
| miniCBF / CBF | `DIRECTORY` + `TEMPLATE` with `#` placeholders, e.g. `sweep_####.cbf`. |
| SMV (ADSC) | Same `DIRECTORY` / `TEMPLATE` / `--image` pattern. Read as MOSFLM does: `XDET` runs along the file's slow axis (DIALS `XDET` is the fast axis). |

For multi-frame CBF/SMV (`FINDSPOTS … MULTIFRAME`), each frame is read from
the `--image` file's sibling with that frame number (same zero padding), else
from `TEMPLATE`; it is an error if neither exists. HDF5 masters are fine as a
single `--image`. A reflection seen on several frames becomes one spot at its
intensity-weighted position and φ (`MERGE_XY 0` keeps one spot per frame;
see `docs/XDS_BRIDGE.md`).

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

Without `DMAX` (or a two-value `RESOLUTION`), `AUTOINDEX` sets the largest
cell edge it searches for as DIALS `find_max_cell` does. With a `CELL`, it is
1.3 × the longest edge. Otherwise it is estimated from the spots:

1. Repeats of a reflection on adjacent frames (within 3 px and 0.5°) are
   merged into one centroid.
2. Reflections at hexagonal-ice ring resolutions (±0.002 Å⁻² in 1/d²) and
   overlapping reflections are dropped.
3. Within 45° blocks of φ, and separately on each side of the rotation axis,
   each reciprocal-lattice point's nearest neighbour gives a direct-space
   distance.
4. The distances are histogrammed (about 5 per bin, longest 1% dropped), and
   1.3 × the far edge of the last bin above 25% of the tallest is the estimate.

Either way it only raises the 200 Å default. It is
logged as `AUTOINDEX: nearest-neighbour max cell … → DMAX …`.

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
| `CELLFIX ON\|OFF\|HARD\|SOFT [s]` | Project the cell onto the space-group metric before PROCESS (default ON, hard). The orientation matrix is rebuilt for the projected cell, keeping its orientation, so predictions use the cell that is reported. |
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
| `REFINE ALL` or `REFINE ORIENTATION CELL DISTANCE BEAM TILT` (any subset) | Parameters positional refinement may change during PROCESS. Default `ALL`; `--xparm` sets `ORIENTATION` (as XDS INTEGRATE; the imported cell and detector are kept). See [Geometry refinement during PROCESS](#geometry-refinement-during-process) for which are refined when. |
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
`CHARACTER n` selects the lattice solution with that lattice character (1–44,
the `character` of each row in `process.json` `indexing.solutions`), and
`SYSTEM tP` selects the lowest-penalty solution of that lattice.
Indexing details live in [`processing-lab/docs/XDS_BRIDGE.md`](processing-lab/docs/XDS_BRIDGE.md).

### Strategy (`RunStep.STRATEGY`)

After indexing (or `XPARM` with the orientation), `STRATEGY` works out the
smallest rotation that records a target completeness, as MOSFLM's
`STRATEGY AUTO` does. Every reflection to the resolution limit is followed
over a full turn: both Ewald crossings are solved, kept when they land on the
detector, and reduced to the asymmetric unit with their Friedel sign. For
each total rotation the start angle recording the most is found.

| Keyword | Meaning |
| ------- | ------- |
| `STRATEGY [AUTO] COMPLETENESS pct` | Completeness wanted, % (default 99). The rotation reported is the smallest reaching it. |
| `STRATEGY … ANOMALOUS` / `NOANOMALOUS` | Aim for anomalous completeness: acentric pairs with both Friedel mates. |
| `STRATEGY … FRAMES first last` | Choose within these collected frames, and report the frames to integrate. `FULL` (default): plan on a full turn. |
| `STRATEGY … ROTATE r` | Also report the best start for a rotation of *r*° (MOSFLM `ROTATE`; `ANGLE` and `RANGE` are aliases). |
| `STRATEGY … RESOLUTION d` | High-resolution limit (Å). Default: `DMIN` if given, otherwise the nearest detector edge. Beyond the edge only the corners record data, so no rotation could be complete. |
| `STRATEGY … STEP deg` | φ bin and rotation step (default 1°). |

The space group matters: higher symmetry needs less rotation. On the CCP4
insulin frame, 95% completeness to 1.87 Å needs 13° in P4₁32, 16° in P23
and 151° in P1, and 95% of anomalous pairs needs 17°, 28° and 218°. The log
prints the best completeness for each rotation, then the chosen segment's
start, end, number of images (and frames) and its statistics by resolution
shell. `process.json` has them under `strategy` (see
[RUN_RECORDS.md](processing-lab/docs/RUN_RECORDS.md)).

```text
XPARM 02-index/XPARM.XDS
SYMMETRY P41212
STRATEGY COMPLETENESS 99 FRAMES 1 1800
GO
```

### Symmetry correction (`RunStep.CORRECT`)

After `PROCESS` has integrated a sweep, `CORRECT` screens axial intensities
for screw-axis periodicity and, for tetragonal data, checks the additional
twofold operators against measured intensities. It recognizes the
P 41 21 2 / P 43 21 2 axial signature, applies P 41 21 2 by convention, and
writes a reindexed `correct.mtz`, leaving `process.mtz` as the integration
output. `correct.log` repeats PROCESS's resolution-bin table (completeness,
R-factors, CC(1/2)) with observations merged and absences removed in the
corrected space group. The opposite enantiomorph cannot be distinguished
from these data alone. Other crystal systems currently get axis scores but no
automatic space-group promotion. The correction is a hypothesis for later
validation; it does not repeat integration.

```text
PROCESS 1 TO 900
GO
CORRECT
GO
```

### Integration

| Keyword | Meaning |
| ------- | ------- |
| `SEPARATION x [y]` | Measurement-box size (mm). |
| `PEAKHALF n` | Peak half-width (pixels); clamped to half `SEPARATION`. Without it the box is sized from the strong spots at the start of PROCESS. |
| `APERTURE ON\|OFF` | Aperture correction (default ON): intensities are divided by the fraction of a spot's flux the measured peak box holds, measured on strong spots by distance from the rotation axis (spot tails along the rotation direction lengthen with it). The log's `APERTURE CORRECTION` line gives the fractions; 1.0 for compact spots (CCP4 insulin), about 0.5–0.7 on 76c1b. Not applied with `PEAKHALF`. |
| `MINFRAC` | Drop reflections with `FRACTIONCALC` below this (default 0.5). |
| `PARTIALSCALE ON\|OFF` | Scale `I`/`σ` by `1/FRACTIONCALC` (default ON). |
| `POLARIZATION ρ` | LP polarisation fraction (default 0.98). |
| `AIR c` | Air absorption (mm⁻¹). |
| `AUTOSIGMA ON\|OFF` | Estimate mosaic / divergence from partials (default ON). |
| `POSTREFINE [PERIOD n] \| OFF` | Misset refine every *n* frames (default period 25). |
| `WIDEN_EWALD ON\|OFF` | Widen Ewald with mosaic (PROCESS default OFF). |
| `NMASKS nx ny` | Standard-profile detector bins. |

### Geometry refinement during PROCESS

Positional refinement fits the predicted spot positions, and where observed
the φ centroids, to strong reflections. Which parameters move depends on how
much rotation the observations cover:

- **Pre-refinement** (two cycles on the first block, before any frame is
  kept): everything `REFINE` allows.
- **Each block** (`POSTREFINE PERIOD`, 25 frames): the misset, the beam and
  the distance only, as XDS INTEGRATE does by default. Over a few degrees the
  spot radii fix only distance ÷ cell, and tilt/twist trades with the misset,
  so a free cell and detector wandered together from block to block (on
  76c1b distance 385 → 392 mm, c 71.8 → 74.5 Å) and pulled the predicted φ
  off the spots. If a block's step exceeds the safety caps, it is retried
  with the misset and beam alone, so the crystal is still followed.
- **The sweep so far** (every 4 blocks, once 15° are covered): the cell, the
  misset and tilt/twist against up to 400 observations from every block,
  pooled — like XDS CORRECT, from a rotation range wide enough to fix the
  cell's shape.

Each `POSITIONAL REFINE` line gives the residual split into detector
(`xy`, px, with its mean offset) and φ (degrees).

### Untrusted pixels (DEFPIX)

PROCESS runs DEFPIX once, before integrating: it samples frames spread over
the range and classifies every pixel, as XDS's DEFPIX does with INIT's
background. Ignored pixels are left out of the background, and a
reflection with one under its peak is rejected.

| Keyword | Meaning |
| ------- | ------- |
| `UNTRUSTED RECTANGLE x1 x2 y1 y2` | Ignore a rectangle (XDS pixels: 1-based, x along fast, bounds included). |
| `UNTRUSTED ELLIPSE x1 x2 y1 y2` | Ignore the ellipse in that bounding box. |
| `UNTRUSTED QUADRILATERAL x1 y1 … x4 y4` | Ignore a quadrilateral (corners in order). |
| `UNTRUSTED CLEAR` | Forget the regions given so far. |
| `EXCLUDE RESOLUTION d1 d2` | Ignore a resolution range (Å, either order). |
| `EXCLUDE ICE` | Ignore the usual ice-ring ranges (3.93–3.87 Å, …). |
| `EXCLUDE CLEAR` | Forget the ranges given so far. |
| `DEFPIX ON\|OFF` | Run DEFPIX (default ON). |
| `DEFPIX FRAMES n` | Frames sampled (default 10, at most 50). Shaded, defective and hot pixels need 3 or more. |
| `DEFPIX VALUE_RANGE lo hi` | A block's background relative to its resolution ring, ×10000 (default 7000 30000, as XDS). |
| `DEFPIX SHADOWS ON\|OFF` | Look for shaded regions (default ON). |
| `DEFPIX DEFECTS ON\|OFF` | Look for defective and hot pixels (default ON). |
| `DEFPIX HALO ON\|OFF` | Look for the beam halo (default ON). |
| `DEFPIX FILE name\|OFF` | Control image (default `defpix.cbf`). |

Several `DEFPIX` options can share a line (`DEFPIX FRAMES 20 SHADOWS OFF`).
DEFPIX ignores:

- **Untrusted regions** and **excluded resolution** as given.
- **Shaded** pixels: the background of each block of pixels (8×8 or larger,
  so a block holds enough counts) is compared with the median block of its
  resolution ring, and blocks outside `VALUE_RANGE` are ignored. A beamstop
  that fills whole rings is compared with the rings just outside it.
- **Defective** pixels: overloaded in at least half the sampled frames
  (stuck), or in two or more and near background in the rest (flickering).
- **Hot** pixels: in every sampled frame at least 10 counts and 20× the
  mean background of their ring.
- The **beam halo**: where the beam is not fully blocked it flares around
  the stop, often filling whole rings, so the shaded test (which compares a
  block with its own ring) cannot see it. DEFPIX walks out from the beam
  along the radial background (each ring's mean without its brightest 5%,
  rings widened until they hold enough counts) from its peak to where the
  excess over the outer background has fallen to a tenth, and ignores every
  pixel inside that radius. A halo must peak at least 1.5× the outer level.
  Integration's low-resolution limit then sits 10 px outside it (e.g.
  `Low-resolution limit for integration: 40.61 Å`), so no measurement box
  reaches the remaining slope, unless `DMAX` or `RESOLUTION` sets the limit.
  The indexing DMAX is unchanged.

FINDSPOTS (over three or more frames) first runs the shaded and halo tests
on five of its frames and leaves those pixels, and any untrusted region, out
of the search: out of the background statistics as well as the peaks, so
the flare neither becomes spots nor raises the threshold.

Defective and hot pixels are only looked for when the frames span at least
5°: over a shorter range a low-resolution spot can stay bright, or
saturated, in every sampled frame.

The background is each pixel's mean over the sampled frames without its
largest value, so a spot crossing a pixel in one sampled frame doesn't count.
The log's `DEFPIX` block lists the counts and the first defective and hot
pixels; `process.json` has them under `integration.defpix`. `defpix.cbf`
holds the background (summed over the sampled frames) for trusted pixels,
−3 for ignored ones and −1 for pixels the detector masks, like XDS's
`BKGPIX.cbf`.

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
| `UNTRUSTED_RECTANGLE` / `_ELLIPSE` / `_QUADRILATERAL` | `UNTRUSTED RECTANGLE` / … |
| `EXCLUDE_RESOLUTION_RANGE` | `EXCLUDE RESOLUTION` |
| `VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS` | `DEFPIX VALUE_RANGE` |

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

`--xparm` (script mode or `run`; or CONTROL `XPARM` / `GXPARM`) seeds:

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
| `process.json` | created at start, updated at every `GO` | `status` (`running` / `completed` / `failed`), `stages` with status and timing, `indexing.selected_solution`, `refinement` (after PROCESS), `integration`, `quality` (unscaled shells + total, and counts that add up to the integrated reflections), `validation`, `warnings`, `error`, `outputs` (relative to `run_directory`). Written atomically; schema version 3 (stages named by `RunStep`: `find_spots`, `indexing`, `defpix`, `integrate`, `strategy`), defined in [processing-lab/docs/RUN_RECORDS.md](processing-lab/docs/RUN_RECORDS.md). |
| `process.md` | with `process.json` | Report rendered from the two records. |
| `process.log` | always (unless `--no-log`) | Everything printed to the terminal. |
| `HKLOUT` (`process.mtz` if unset) | PROCESS `GO` | Unmerged MOSFLM-style MTZ. |
| `spots.spt` (or `FINDSPOTS FILE`) | FINDSPOTS `GO` | MOSFLM spot list: header (nfast, nslow, pixel, yscale, omega; invflag; beam mm) then `x_mm y_mm φ_start φ_end I σ`. |
| `XPARM.XDS` | successful AUTOINDEX | Bridge-format orientation + cell + detector (written next to `DIRECTORY`). |
| `GXPARM.XDS` | PROCESS `GO` (unless `XPARMOUT OFF`) | Same layout after integration. Override the path with `XPARMOUT`. |
| `defpix.cbf` (or `DEFPIX FILE`) | PROCESS `GO` (unless `DEFPIX OFF` / `DEFPIX FILE OFF`) | [DEFPIX](#untrusted-pixels-defpix) control image: background for trusted pixels, −3 ignored, −1 masked by the detector. |
| `index.log` | successful AUTOINDEX | Autoindexing result: `CELL` / `MATRIX` / `SYSTEM` / `CHARACTER` / `NSPOTS` / `INDEXED` / `RMSD_HKL`. |
| `xds/`, `xia2/`, `mosflm/` | start (with an output directory); `xds/` gets the XPARMs after the run | [Repeat inputs](#repeat-inputs-xds-xia2-mosflm) for each run mode. |
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
