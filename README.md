# ProcessingLab

Progressive C++ port of the MOSFLM diffraction-processing core.

![ProcessingLab Logo](resources/ProcessingLabLogo.png)



## Build and install

Requirements: CMake ≥ 3.24, a C++20 compiler, ZLIB and a CCP4 installation
(`configure` reads `$CCP4`). Optional: gfortran, for DPS autoindexing; HDF5,
for Eiger/NeXus masters; and FFTW.

```bash
./configure --prefix=$HOME/.local   # or --with-ccp4-prefix=/Applications/ccp4-9
make                                # → bin/ProcessingLab
make test                           # ctest in processing-lab/build
make install                        # → $(prefix)/bin/ProcessingLab
```

- **`./configure`:** without `--prefix`, it installs into the CCP4 tree
  (`$CCP4/bin`). `--with-ccp4-prefix` names the CCP4 tree when `$CCP4` is
  not set.
- **`make`:** builds CBFlib, the DPS archives (skipped with a warning when
  gfortran is missing; then `make cpp-core` builds without autoindexing) and
  `ProcessingLab`.
- **`make install`:** copies only what `make` built, and stops if you
  haven't run `make` yet. Run `make` as yourself; if the prefix needs root,
  use `sudo` only for `make install`, so nothing in the tree ends up owned
  by root.
- **Removing it:** `make uninstall` takes it out again.

The developer tool `processinglab_dev` stays in `processing-lab/build/` and
is never installed. If you configure `processing-lab` with CMake directly,
build CBFlib first: `make -C cbf`.

```bash
./bin/ProcessingLab --version
./bin/ProcessingLab HKLOUT /tmp/out.mtz --image /path/to/insulin_1_001.img \
  processing-lab/testdata/images/known_indexed/cpp_process.in
# C++ vs F77 PROCESS wall-time (CCP4 ipmosflm by default):
./processing-lab/tools/bench_ipmosflm_perf.sh
```

Top-level CMake (`-S .`) adds `processing-lab` when `PROCESSINGLAB_BUILD_CPP=ON`
(default).

## Optional Fortran oracle

Requires CCP4. See [`porting.md`](porting.md) for macOS notes.

```bash
make fortran-oracle   # → bin/ipmosflm (F77; does not overwrite ProcessingLab)
# or:
#   ./configure --with-ccp4-prefix=/Applications/ccp4-9
#   make -C mosflm MOSFLM_FORTRAN_ORACLE=1 ipmosflm
# Experimental C++ PREDICT splice into the oracle:
#   cmake -S processing-lab -B processing-lab/build -DMOSFLM_FORTRAN_ORACLE=ON
#   cmake --build processing-lab/build --target processinglab_predict_hook
#   make -C mosflm MOSFLM_FORTRAN_ORACLE=1 MOSFLM_CPP_PREDICT=1 ipmosflm
```

# Usage

There are two ways to run ProcessingLab:

- **`ProcessingLab run`**: give it the images and it finds spots, autoindexes
  and integrates, taking the experiment metadata from the image headers.
- **CONTROL script**: MOSFLM-style keywords on stdin or in a file, for full
  control of every step.

`run` always writes into a new output directory, and a script does too with
`--output-dir`. The full keyword and option reference is in
[`USAGE.md`](USAGE.md).

## `ProcessingLab run`

```bash
# Whole sweep from an Eiger HDF5 master
./bin/ProcessingLab run --image ~/Downloads/thaumatin_9_1_master.h5

# CBF sweep: name any one frame; its numbered siblings are the sweep
./bin/ProcessingLab run --image ../test_xray_data_cbf/76c1b_0002212_1_0001.cbf \
    --range 1:900 --output-dir results/76c1b

# Take the cell, symmetry and geometry from XDS.INP (or .xinfo)
./bin/ProcessingLab run --image ../test_xray_data_cbf/76c1b_0002212_1_0001.cbf \
    --xds-inp XDS.INP --range 1:900

# An .xinfo with several sweeps needs --sweep
./bin/ProcessingLab run --xinfo automatic.xinfo --sweep SWEEP1

# State where the parameters come from: here, only XDS.INP
./bin/ProcessingLab run --run-mode xds --xds-inp XDS.INP

# Re-process with the refined cell, geometry and orientation of a run
./bin/ProcessingLab run --from processing_20261022_235027
```

| Option | Role |
| ------ | ---- |
| `--image PATH` | HDF5/NeXus master, or one frame of a CBF/SMV sweep. Repeatable: each image is one sweep, processed in turn into `sweep_1/`, `sweep_2/`, … |
| `--range FIRST:LAST` | Inclusive frame range. Default: every frame found. |
| `--output-dir DIR` | Output directory; it must not exist yet. Default: a new `processing_YYYYMMDD_HHMMSS`. |
| `--xds-inp PATH` | Import `XDS.INP`. |
| `--xinfo PATH` | Import a Xia2 `.xinfo`. |
| `--sweep NAME` | Which `.xinfo` sweep, when it has several. |
| `--xparm PATH` | Import `XPARM.XDS` / `GXPARM.XDS`: cell, space group, geometry and orientation. With an orientation, autoindexing is skipped. |
| `--from RUN` | Re-process from a completed run (its directory or `process.json`), with its refined parameters; see below. |
| `--spacegroup N` | Space group number or symbol, over every import and `--from`. |
| `--cell A B C α β γ` | Unit cell, over every import and `--from`; must fit the space group. |
| `--run-mode MODE` | Where the parameters come from, recorded in the log, `input.json` and `process.md`: `processinglab` (`RunMode.PROCESSINGLAB`, the default: image headers, overridden by any imports), `xds` (`RunMode.XDS`: needs `--xds-inp`, no `--xinfo`), `xia2` (`RunMode.XIA2`: needs `--xinfo`, no `--xds-inp`). `mosflm` (`RunMode.MOSFLM`) is a CONTROL script run (see below). |

### Run modes (`--run-mode`)

`--run-mode` states where a run's parameters come from, so a log or report
shows it at a glance. ProcessingLab refuses inputs that contradict the mode.

| Mode | Shown as | Parameters from | Requires |
| ---- | -------- | --------------- | -------- |
| `processinglab` | `RunMode.PROCESSINGLAB` | Image headers, overridden per field by any `--xds-inp` / `--xinfo` / `--xparm` / `--from` | Nothing (default for `run`) |
| `xds` | `RunMode.XDS` | `XDS.INP`; headers fill unset fields | `--xds-inp`; no `--xinfo` / `--sweep` / `--xparm` / `--from` |
| `xia2` | `RunMode.XIA2` | Xia2 `.xinfo`; headers fill unset fields | `--xinfo`; no `--xds-inp` / `--xparm` / `--from` |
| `mosflm` | `RunMode.MOSFLM` | MOSFLM CONTROL keywords | A CONTROL script, not `run` (default for scripts) |

```bash
./bin/ProcessingLab run --image sweep_0001.cbf                    # RunMode.PROCESSINGLAB
./bin/ProcessingLab run --run-mode xds --xds-inp XDS.INP          # RunMode.XDS
./bin/ProcessingLab run --run-mode xia2 --xinfo automatic.xinfo   # RunMode.XIA2
./bin/ProcessingLab --run-mode mosflm --image img.cbf script.in   # RunMode.MOSFLM
```

The mode is recorded in `process.log` (`Run mode:  RunMode.XDS (parameters
from …)`), in `input.json` (`run_mode`, `parameter_source`) and in
`process.md`. A mismatch fails before any processing, with what to change:

```text
ProcessingLab run: RunMode.XDS needs --xds-inp XDS.INP
ProcessingLab: RunMode.XDS is for 'ProcessingLab run'; a CONTROL script is RunMode.MOSFLM
```

### Metadata precedence

Metadata precedence, per field, lowest to highest: image header, `XDS.INP`,
`.xinfo`, XPARM (`--xparm`, or `--from`'s refined values), command line. A
source only overrides the fields it sets.
`process.log` opens with each value and where it came from, and warns when
sources disagree:

```text
 Experiment metadata
  Image source:  .../76c1b_0002212_1_0001.cbf  [CLI]
  Wavelength:    0.97949 Å  [XDS.INP]
  Distance:      384.99 mm  [XDS.INP]
  Beam:          214.28 207.96 mm  [image header]
  Cell:          62.04 186.93 72.36 90.00 100.19 90.00  [XDS.INP]
  Frame range:   1–900  [CLI]  (900 on disk)
```

Spots for indexing come from the first 50 frames of the range. If you import
a cell, the longest cell edge allowed during autoindexing (`DMAX`) is set to
1.7× its longest edge. Otherwise it stays at the 200 Å default, which is too
short for very long axes (76c1b's 187 Å). For those, pass `--xds-inp` or use a
script with `DMAX`.

### Overriding the space group and cell (`--spacegroup`, `--cell`)

To use a space group or cell you've determined yourself, give them on the
command line. They win over every other source, `--from` included:

```bash
# Thaumatin is P4₁2₁2 (92), not P4₁ (76): keep the refined cell, change the symmetry
./bin/ProcessingLab run --from processing_20261022_235027 --spacegroup 92

# A cell of your own as well (six numbers, or one quoted list)
./bin/ProcessingLab run --from processing_20261022_235027 --spacegroup P41212 \
    --cell 57.898 57.898 150.952 90 90 90
```

- **`--spacegroup`:** takes a number or a symbol (`92`, `P41212`).
- **`--cell`:** must fit the space group's crystal system (tetragonal:
  a = b, α = β = γ = 90°). A cell that doesn't fit is refused before
  anything runs.
- **With `--from` or `--xparm`:** the orientation matrix also encodes the
  cell, so it's rebuilt with your cell and the refined orientation kept;
  otherwise predictions would keep the old cell. The cell then stays fixed
  during integration, as the orientation alone is refined.
- **In the records:** each override is tagged `[CLI]` in `process.log` and
  `input.json`, with a warning naming the value it replaced. It also goes
  into the repeat inputs (`xds/`, `xia2/`, `mosflm/`).

Changing only the space group within the same lattice (P4₁ → P4₁2₁2) leaves
the predictions untouched. On thaumatin frames 1–100 the residual stays at
0.44–0.55 px. A cell from another source, such as a PDB `CRYST1`, may not
match this dataset at the refined distance. With 57.898, 57.898, 150.952 Å
against the refined 57.938, 57.938, 150.244 Å, the residual rose from 0.47
to 1.12 px. Prefer `--spacegroup` alone, and keep a reference cell for
scaling or refinement unless the data agree with it.

### Repeat inputs (`xds/`, `xia2/`, `mosflm/`)

At the start of every run, ProcessingLab writes the resolved experiment as
input for each run mode, so the run can be repeated, or handed to XDS or
xia2:

| Directory | Files | Repeat with |
| --------- | ----- | ----------- |
| `xds/` | `XDS.INP`; after the run, copies of `XPARM.XDS` and `GXPARM.XDS` | `ProcessingLab run --run-mode xds --xds-inp RUN/xds/XDS.INP` |
| `xia2/` | `automatic.xinfo`, `run_xia2.sh` | `ProcessingLab run --run-mode xia2 --xinfo RUN/xia2/automatic.xinfo` |
| `mosflm/` | `mosflm.in`: the run's CONTROL script, without `HKLOUT` | `ProcessingLab --run-mode mosflm --output-dir NEW RUN/mosflm/mosflm.in` |

Each reproduces the original run. Insulin (1 frame) and thaumatin
(20 frames) give the same reflections and statistics in all three modes.
CONTROL script runs with `--output-dir` write the same three directories.

- **`XDS.INP`:** for XDS itself, it is a complete starting point:
  - `JOB`, `DETECTOR` and `OVERLOAD` come from the header;
  - the spot and background ranges match the run;
  - an Eiger sweep appears as `x_??????.h5`, and needs XDS's HDF5 `LIB=`
    plugin;
  - check `DETECTOR`, `ROTATION_AXIS` and the trusted region for your
    beamline.
- **The copied `XPARM.XDS` and `GXPARM.XDS`:** these use ProcessingLab's
  own layout. Use them with `--xparm`, not with XDS.
- **`automatic.xinfo`:** it holds `USER_CELL` and `USER_SPACEGROUP` when the
  run knew them. `BEAM` and `DISTANCE` appear only where the run overrode
  the image header; xia2 reads the header for the rest.
- **Refined values:** these files hold the values at the start of the run.
  To re-process from the refined values, use `--from` (above).

### Re-processing with refined parameters (`--from`)

After a run has refined the cell and geometry, process again starting from
them with `--from`. It takes a finished run directory, or its `process.json`:

```bash
./bin/ProcessingLab run --image ~/Downloads/thaumatin_9_1_master.h5 --range 1:900
./bin/ProcessingLab run --from processing_20261022_235027 --output-dir thaumatin_rerun
```

- **What `--from` takes:** the refined cell, space group, wavelength,
  distance, beam, pixel size and orientation from that run's `GXPARM.XDS`,
  and the refined detector tilt/twist from its `process.json`. `GXPARM.XDS`
  doesn't carry tilt/twist. Without them, thaumatin's positional residual
  rose from 0.5 to 2 px.
- **Images and frame range:** taken from its `input.json`, unless
  `--image` or `--range` names others.
- **No autoindexing:** the orientation is known, so `AUTOINDEX` is skipped
  and only the orientation is refined. The cell, distance and beam stay at
  their refined values.
- **Provenance:** each value is tagged `[previous run]` in `process.log`,
  and `input.json` records the run in `inputs.imports.previous_run`.
- **What it refuses:** a run that failed, is still running or never reached
  `PROCESS`. For a multi-sweep run, give one `sweep_N` directory.

`--xparm GXPARM.XDS` does the same from an XPARM file alone (for example,
one from XDS). It leaves out the tilt/twist that only `process.json` holds.

## CONTROL script

```bash
./bin/ProcessingLab --image ~/Downloads/thaumatin_9_1_master.h5 \
    --output-dir thaumatin_run <<'EOF'
FINDSPOTS THRESHOLD 5.0 MAXSPOTS 1500 FIND 1 FILE thaum.spt MULTIFRAME RANGE 1 50
GO
AUTOINDEX DPS THRESHOLD 15
GO
PROCESS 1 TO 999
GO
EXIT
EOF
```

```bash
./bin/ProcessingLab --image ../test_xray_data_cbf/76c1b_0002212_1_0001.cbf \
    --output-dir 76c1b_run <<'EOF'
DIRECTORY ../test_xray_data_cbf
TEMPLATE 76c1b_0002212_1_####.cbf
DMAX 320
FINDSPOTS THRESHOLD 5.0 MAXSPOTS 1500 FIND 1 FILE spots.spt MULTIFRAME RANGE 1 50
GO
AUTOINDEX DPS THRESHOLD 15 SEPMIN 6 CLUSTER_RADIUS 2
GO
PROCESS 1 TO 900
GO
EXIT
EOF
```

- **Script source:** CONTROL is read from a file named on the command line,
  otherwise from stdin; `GO` runs the pending step and `EXIT` ends the session.
- **`--output-dir DIR`:** the directory must be new. Relative paths in the
  script (`DIRECTORY`, `IMAGE`) still resolve against the current directory.
  Without `--output-dir`, the log, MTZ and `ProcessingLabIndex.txt` go in
  the current directory, and spot lists, `XPARM.XDS` and `GXPARM.XDS` go in
  `DIRECTORY`.
- **Run mode:** a script run is `RunMode.MOSFLM`: parameters come from its
  MOSFLM keywords. `--run-mode mosflm` says so explicitly; other modes are
  refused here.
- **Imports:** `--xds-inp`, `--xinfo` and `--xparm XPARM.XDS` fill CONTROL
  before the script runs, in that order, and the script's keywords override
  them.
- **Threads:** by default all cores are used. `NPROC n` sets the number of
  threads; the log reports it as `USING n PROCESSORS`.
- **CBF/SMV sweeps:** give `DIRECTORY` and `TEMPLATE` (`#` for the frame
  digits). An HDF5 master needs only `--image`.

## Outputs

| File | Contents |
| ---- | -------- |
| `input.json` | What ProcessingLab was asked to process, resolved before it starts: images, imports, every experiment value with its source and what it overrode, and the effective settings |
| `process.json` | What it found: stage status, the selected indexing solution, refined geometry, integration and resolution-shell statistics, warnings, errors and output paths |
| `process.md` | A readable report generated from the two JSON files |
| `process.log` | Everything printed, starting with the CONTROL input (and, for `run`, the resolved metadata) |
| `process.mtz` | Unmerged MTZ, for POINTLESS/AIMLESS (`HKLOUT` renames it) |
| `spots.spt` (or the `FILE` name) | Spot list from `FINDSPOTS` |
| `XPARM.XDS` | Orientation, cell and detector after `AUTOINDEX` |
| `GXPARM.XDS` | The same after `PROCESS` (refined) |
| `ProcessingLabIndex.txt` | Autoindexing result in a simple text form |
| `xds/`, `xia2/`, `mosflm/` | Inputs to repeat the run in each run mode (see [Repeat inputs](#repeat-inputs-xds-xia2-mosflm)) |

`run` always writes the three run records; a CONTROL script writes them when
given `--output-dir`. Each value lives in one place: inputs and their sources
in `input.json`, results in `process.json` (`process.md` is generated from
them), per-image detail in `process.log`, reflections in `process.mtz`.

`process.json` is created as `"status": "running"` with every stage
`pending`. Each stage then goes `running` → `completed` or `failed`, and the
run ends `completed` or `failed` with the error. Every update replaces the file
atomically, so a reader never sees a partial file. A run that crashed is left
at `running`, showing the stage it was in. Both files carry
`"schema_version": 1`.

## Test fixtures

Some differential tests compare against `ipmosflm` runs on the CCP4 insulin
example. The fixture MTZs are not in git. Regenerate them with
`processing-lab/tools/capture_ipmosflm_oracle.sh`, which sets the macOS
library path itself and prints the end of the `ipmosflm` log if a step fails:

- `testdata/oracle/known_indexed/insulin_1_001.mtz`: run the script with
  CCP4's `ipmosflm` (the default when `bin/ipmosflm` does not exist).
- `testdata/oracle/known_indexed_cpp_predict/insulin_1_001.mtz`: build
  `ipmosflm` with `MOSFLM_CPP_PREDICT=1` (see above), run the script with
  `IPMOSFLM=` pointing at it, and copy the MTZ there.

# Binaries
**`ProcessingLab`** (`ProcessingLab.exe` on Windows) — a thin CONTROL →
PROCESS driver that does **not** link the Fortran amalgamation. Autotools
`mosflm/` remains an optional scientific oracle named **`ipmosflm`** (F77).

| Path | Role |
| ---- | ---- |
| [`processing-lab`](processing-lab/) | C++ geometry, PREDICT, integration, CONTROL, CTest |
| [`mosflm/`](mosflm/) | Amalgamated Fortran (**opt-in oracle** → `bin/ipmosflm`) |
| `bin/ProcessingLab` | C++ CONTROL binary (default; `.exe` on Windows) |
| `processing-lab/build/processinglab_dev` | Developer diagnostics (not installed) |
| `bin/ipmosflm` | F77 oracle only (`make fortran-oracle`) |

## License

BSD-style — see [`COPYING`](COPYING) (A.G.W. Leslie / MRC LMB).
