# ProcessingLab

**A C++20 diffraction-processing engine for macromolecular crystallography.**

![ProcessingLab Logo](resources/ProcessingLabLogo.png)

ProcessingLab is a modern, progressively independent C++ implementation of
the MOSFLM diffraction-processing core. It finds spots, autoindexes, refines
the orientation and integrates diffraction images. It takes its experiment
parameters from image headers, `XDS.INP`, Xia2 `.xinfo` files and
MOSFLM-style CONTROL scripts.

ProcessingLab is a component of **DiffractionLab**, but its processing engine
and command-line interface can also be used on their own.

## Features

- **Automated processing:** give it the images, and it resolves the
  experiment metadata, finds spots, autoindexes and integrates the
  reflections.
- **MOSFLM-compatible CONTROL interface:** explicit control of each
  processing stage, with familiar keywords and `GO`.
- **Image formats:** CBF and SMV sweeps, and HDF5/NeXus master files
  (Dectris Eiger).
- **Metadata import and provenance:** experiment parameters from image
  headers, `XDS.INP`, Xia2 `.xinfo` and XDS orientation files. Every value
  records where it came from.
- **Crystallographic constraints:** unit cells are checked against the
  space group, and supplied constraints are kept separate from refined
  results.
- **Reproducible processing:** each run records its effective inputs and
  their sources, its stages, warnings, results and outputs.
- **Interoperability:** unmerged MTZ, XDS-style orientation files, and
  generated inputs to repeat a run with XDS, Xia2 or a MOSFLM-style script.
- **C++ processing core:** geometry, reflection prediction, integration and
  CONTROL execution are native C++20. DPS autoindexing optionally uses
  Fortran.
- **Differential testing:** selected results are compared against the
  original Fortran `ipmosflm`.

ProcessingLab covers indexing and integration. **Scaling and merging happen
downstream**, typically with CCP4 Pointless and Aimless.

## Installation

### Requirements

| Dependency | Purpose |
| ---------- | ------- |
| CMake ≥ 3.24 | Build system |
| C++20 compiler | Processing engine |
| ZLIB | Compression |
| CCP4 | Required by `configure` (`$CCP4`, or `--with-ccp4-prefix`); default install location; reference `ipmosflm` |
| gfortran (optional) | DPS autoindexing |
| HDF5 (optional) | Eiger / NeXus master files |
| FFTW (optional) | FFT-based autoindexing |

### Build and install

```bash
./configure --prefix="$HOME/.local"   # or --with-ccp4-prefix=/Applications/ccp4-9
make                                  # → bin/ProcessingLab
make test                             # builds and runs the CTest suite
make install                          # → $prefix/bin/ProcessingLab
```

- **The prefix:** without `--prefix`, ProcessingLab installs into the CCP4
  tree (`$CCP4/bin`). `--with-ccp4-prefix` names the CCP4 installation when
  `$CCP4` is not set.
- **Without gfortran:** `make` skips the DPS archives with a warning, and
  `make cpp-core` builds the C++ core alone.
- **Permissions:** run `make` as yourself. If the prefix needs elevated
  permissions, use `sudo` for `make install` only: it copies what `make`
  built and never compiles.
- **Removing it:** `make uninstall`.

Check the installation:

```bash
./bin/ProcessingLab --version
```

Direct CMake builds, the optional Fortran reference implementation and the
test fixtures are described in [DEVELOPMENT.md](DEVELOPMENT.md).

## Quick start

### 1. Automated processing

An HDF5/NeXus master file:

```bash
./bin/ProcessingLab run --image ~/data/thaumatin_9_1_master.h5
```

ProcessingLab reads the experiment metadata from the image headers, finds
spots, autoindexes and integrates, all in a new output directory.

A CBF sweep: name any one frame, and optionally an inclusive frame range:

```bash
./bin/ProcessingLab run \
    --image ../test_xray_data_cbf/76c1b_0002212_1_0001.cbf \
    --range 1:900 \
    --output-dir results/76c1b
```

With an experiment description from XDS:

```bash
./bin/ProcessingLab run \
    --image ../test_xray_data_cbf/76c1b_0002212_1_0001.cbf \
    --xds-inp XDS.INP \
    --range 1:900
```

Or from Xia2. An `.xinfo` with several sweeps needs `--sweep`:

```bash
./bin/ProcessingLab run --xinfo automatic.xinfo --sweep SWEEP1
```

### 2. Reprocess a completed run

Start from the refined parameters of an earlier run:

```bash
./bin/ProcessingLab run \
    --from processing_20261022_235027 \
    --output-dir thaumatin_rerun
```

`--from` starts from the earlier run's final model: the refined cell, space
group, wavelength, detector geometry (including tilt and twist), mosaicity
and orientation. The model is checked first. A matrix that doesn't belong to
its cell, a cell its space group forbids, or a wavelength that disagrees with
the earlier inputs is refused.

With the orientation known, autoindexing is skipped. By default only the
orientation is refined and the mosaicity re-estimated; `--refine` chooses
otherwise, for example `--refine orientation,cell`. Give a completed run directory or its `process.json`; for a
multi-sweep run, the `sweep_N` directory you want.

### 3. A MOSFLM-style CONTROL script

For explicit control of each stage:

```bash
./bin/ProcessingLab \
    --image ~/data/thaumatin_9_1_master.h5 \
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

- **Script source:** CONTROL comes from a file named on the command line,
  or from standard input.
- **Running:** `GO` runs the pending step, and `EXIT` ends the session.
- **Images:** CBF and SMV sweeps are named with `DIRECTORY` and `TEMPLATE`.
  An HDF5 master needs only `--image`.

[USAGE.md](USAGE.md) is the complete command-line and CONTROL keyword
reference.

## What a run produces

From `run --image thaumatin_9_1_master.h5 --range 1:100`:

```text
processing_20261023_094808/
├── input.json            resolved inputs, with the source of every value
├── process.json          stage status, results, warnings, outputs
├── process.md            readable report generated from the two records
├── process.log           full processing log
├── process.mtz           unmerged reflections
├── spots.spt             FINDSPOTS spot list
├── index.log             autoindexing result (cell, matrix, lattice, indexed %)
├── XPARM.XDS             orientation, cell and detector after AUTOINDEX
├── GXPARM.XDS            the same after PROCESS (refined)
├── xds/                  XDS.INP, XPARM.XDS, GXPARM.XDS
├── xia2/                 automatic.xinfo, run_xia2.sh
└── mosflm/               mosflm.in
```

Excerpt of its `process.md`:

```markdown
**Status:** completed
**Run mode:** `RunMode.PROCESSINGLAB` — parameters from image headers, …

| Parameter | Value | Source | Overrides |
| --- | --- | --- | --- |
| Wavelength (Å) | 0.979499 | image_header | — |
| Distance (mm) | 200 | image_header | — |
| Frames | 1–100 | command_line | 1–2400 (image_header) |

## Indexing
- Lattice: tP, space group 75
- Indexed 97% of 1500 spots, RMSD 0.093 (hkl)

## Refined geometry
- Cell: 57.938 57.938 150.244 Å, 90.00 90.00 90.00°
- Distance: 199.982 mm; beam 166.119 172.925 mm
- Mosaicity: 0.101°

## Integration
- 60881 reflections from 543070 predictions (44 overlapped, 512506 profile-fitted)

## Data quality
P 4; unscaled: no batch scaling or outlier rejection.
| d_min (Å) | Observed | Unique | Completeness | R_meas | I/σ | CC½ |
| **total** | 58405 | 45327 | 26.7% | 17.9% | 3.30 | 97.6 |
```

The statistics are unscaled and in the indexing space group. They show how
the integration went, not the final data quality, which comes from scaling.

## Command-line interface

`run` is the recommended entry point for automated processing.

| Option | Description |
| ------ | ----------- |
| `--image PATH` | HDF5/NeXus master, or one frame of a CBF/SMV sweep. Repeatable: each image is a sweep, processed into `sweep_1/`, `sweep_2/`, … |
| `--range FIRST:LAST` | Inclusive frame range. Default: every frame found. |
| `--output-dir DIR` | Output directory; it must not exist yet. |
| `--xds-inp PATH` | Import experiment parameters from `XDS.INP`. |
| `--xinfo PATH` | Import experiment parameters from a Xia2 `.xinfo`. |
| `--sweep NAME` | The `.xinfo` sweep to use, when it has several. |
| `--xparm PATH` | Import an `XPARM.XDS` / `GXPARM.XDS`: cell, space group, geometry and orientation. |
| `--from RUN` | Reprocess from a completed run, or its `process.json`. |
| `--spacegroup N` | Override the space group, by number or symbol. |
| `--cell A B C α β γ` | Override the unit cell; it must fit the space group. |
| `--refine LIST` | What integration refines: `orientation`, `cell`, `distance`, `beam`, `tilt_twist`, `mosaicity`, `all` or `none`; the rest stay fixed. |
| `--run-mode MODE` | State where the parameters come from (see below). |

Without `--output-dir`, each run creates a new `processing_YYYYMMDD_HHMMSS`
directory. A directory is never reused, so runs stay independent and
earlier results are never overwritten.

## Run modes and metadata provenance

A run mode states where a run's parameters come from, so it can be
reproduced.

| Mode | Parameters from | Requires |
| ---- | --------------- | -------- |
| `processinglab` (default for `run`) | Image headers, overridden per field by any imports | Nothing |
| `xds` | `XDS.INP`; headers fill unset fields | `--xds-inp` |
| `xia2` | Xia2 `.xinfo`; headers fill unset fields | `--xinfo` |
| `mosflm` (default for scripts) | MOSFLM-style CONTROL keywords | A CONTROL script or standard input |

```bash
./bin/ProcessingLab run --image sweep_0001.cbf                    # RunMode.PROCESSINGLAB
./bin/ProcessingLab run --run-mode xds --xds-inp XDS.INP          # RunMode.XDS
./bin/ProcessingLab run --run-mode xia2 --xinfo automatic.xinfo   # RunMode.XIA2
```

The mode is recorded in `process.log`, `input.json` and `process.md`.
Inputs that contradict it are refused before processing starts.

### Parameter precedence

For `run`, each field is resolved separately, lowest priority first:

1. Image header
2. `XDS.INP`
3. Xia2 `.xinfo`
4. `--xparm`, or the refined values of `--from`
5. The command line (`--range`, `--spacegroup`, `--cell`)

A source overrides only the fields it supplies; it never erases a value it
doesn't set. `process.log` opens with each resolved value and its source,
and warns when sources disagree:

```text
 Experiment metadata
  Wavelength:    0.97949 Å  [XDS.INP]
  Distance:      384.99 mm  [XDS.INP]
  Beam:          214.28 207.96 mm  [image header]
  Cell:          62.04 186.93 72.36 90.00 100.19 90.00  [XDS.INP]
  Frame range:   1–900  [CLI]  (900 on disk)
```

In a CONTROL script run, the imports fill the initial configuration and the
script's keywords override them.

### Unit-cell and space-group constraints

```bash
./bin/ProcessingLab run \
    --from processing_20261022_235027 \
    --spacegroup P41212 \
    --cell 57.898 57.898 150.952 90 90 90
```

- **The cell must fit the space group's crystal system.** A tetragonal cell
  needs a = b and α = β = γ = 90°. An incompatible cell is refused before
  processing.
- **The orientation matrix encodes the cell.** When you override the cell
  while reusing an orientation (`--from`, `--xparm`), the matrix is rebuilt
  with your cell and the refined orientation kept. Otherwise predictions
  would still use the old cell.

Changing only the space group within one lattice, such as P4₁ to P4₁2₁2,
leaves the predictions unchanged. A reference cell from elsewhere, such as a
PDB `CRYST1` record, may not match the metric refined from these images. On
thaumatin, imposing 57.898, 57.898, 150.952 Å over the refined 57.938,
57.938, 150.244 Å raised the positional residual from 0.47 to 1.12 px.
Change only the space group when that is all the evidence supports, and
check the positional residuals in `process.log` before accepting an imposed
cell.

### Autoindexing and long cell axes

The autoindexing search is bounded by `DMAX`, which is 200 Å by default.
With an imported cell, `run` sets it to 1.7× the longest cell edge. For
crystals with very long axes (the 187 Å axis of the 76c1b test data),
supply the cell with `--xds-inp`, or set `DMAX` in a CONTROL script.

A search bound doesn't make an indexing solution correct. Check the selected
lattice and cell against the fraction of spots indexed, the positional
residuals and the symmetry evidence.

## Processing records

| File | Purpose |
| ---- | ------- |
| `input.json` | Resolved inputs, with each value's source and what it overrode, the images and the effective settings |
| `process.json` | Stage status, the selected indexing solution, refined parameters, integration statistics, warnings, errors and output paths |
| `process.md` | Readable report generated from `input.json` and `process.json` |
| `process.log` | Detailed log, including the CONTROL input and resolved metadata |
| `process.mtz` | Unmerged reflections for scaling |
| `spots.spt` | Spot list from `FINDSPOTS` (or the `FILE` name given) |
| `index.log` | Autoindexing result: `CELL`, `MATRIX`, `SYSTEM`, `INDEXED`, `RMSD_HKL` |
| `XPARM.XDS` | Orientation, cell and detector after autoindexing |
| `GXPARM.XDS` | Refined orientation, cell and detector after integration |

Every information type has one home:

- inputs and their provenance: `input.json`;
- results and status: `process.json`;
- per-image diagnostics: `process.log`;
- reflections: `process.mtz`.

The JSON records give applications and pipelines a stable interface. Both
carry `"schema_version": 2`; every field, unit and count is defined in
[processing-lab/docs/RUN_RECORDS.md](processing-lab/docs/RUN_RECORDS.md).
`--from` also reads version 1 runs.

`process.json` starts as `"status": "running"` with every stage `pending`.
Each stage then moves to `running` and on to `completed` or `failed`, and
the run ends `completed` or `failed` with the error. Each update replaces
the file atomically, so a reader never sees it half-written. If the process
dies, the file is left at `running`, showing the stage it was in.

`run` always writes these records. A CONTROL script writes them when given
`--output-dir`.

## Interoperability with XDS, Xia2 and MOSFLM

When a run starts, ProcessingLab writes its resolved experiment as input for
each run mode, so the run can be repeated or handed to another program. A
CONTROL script run with `--output-dir` does the same.

| Directory | Contents | Use |
| --------- | -------- | --- |
| `xds/` | `XDS.INP`; after the run, copies of `XPARM.XDS` and `GXPARM.XDS` | Repeat with `--run-mode xds`, or start an XDS run |
| `xia2/` | `automatic.xinfo`, `run_xia2.sh` | Repeat with `--run-mode xia2`, or run Xia2 |
| `mosflm/` | `mosflm.in`: the run's CONTROL script, without `HKLOUT` | Repeat as a CONTROL-script run |

```bash
./bin/ProcessingLab run --run-mode xds \
    --xds-inp processing_20261022_235027/xds/XDS.INP
./bin/ProcessingLab run --run-mode xia2 \
    --xinfo processing_20261022_235027/xia2/automatic.xinfo
./bin/ProcessingLab --run-mode mosflm --output-dir repeated_run \
    processing_20261022_235027/mosflm/mosflm.in
```

Each repeat reproduces the original run. The insulin and thaumatin tests
give the same reflections and statistics in all three modes. The files hold
the parameters from the start of the original run; to use the refined
parameters, use `--from`.

**Notes for XDS:**

- **`XDS.INP`:** it is a complete starting point: detector type, overload,
  image template, data/spot/background ranges and geometry. The detector,
  overload and rotation axis are derived from the image header, so check
  them, and the trusted region, for your beamline.
- **Eiger HDF5:** the sweep appears as `name_??????.h5`, and needs XDS's
  HDF5 plugin (`LIB=`).
- **XPARM files:** ProcessingLab's `XPARM.XDS` and `GXPARM.XDS` use its own
  layout. Use them with `--xparm`; XDS can't read them.
- **`automatic.xinfo`:** `BEAM` and `DISTANCE` appear only where the run
  overrode the image header; Xia2 reads the header for the rest.

## Scaling and merging with CCP4

ProcessingLab writes unmerged reflections. Pointless then assesses the
space group, and Aimless scales and merges:

```bash
pointless HKLIN process.mtz HKLOUT process_pointless.mtz
aimless HKLIN process_pointless.mtz HKLOUT process_aimless.mtz <<'EOF'
GO
EOF
```

[`doc/reindex_scale.com`](doc/reindex_scale.com) runs the same two steps.
This is a minimal example: a given experiment or CCP4 version may need more
keywords. Read Pointless's symmetry assessment and Aimless's merging
statistics; a run that completes doesn't establish the right symmetry or
good data.

## Testing and validation

```bash
make test
```

CTest runs the unit and differential tests. Some compare ProcessingLab
against the Fortran `ipmosflm` on the CCP4 insulin example. Their reference
MTZs aren't in git; [DEVELOPMENT.md](DEVELOPMENT.md) describes how to
regenerate them.

For scientific validation, look beyond whether processing completes. Check:

- the refined cell and lattice;
- the fraction of spots indexed;
- the positional residuals;
- completeness, I/σ and the resolution-dependent statistics.

Evaluate an imposed cell or space group against the diffraction evidence.

## Project layout

| Path | Role |
| ---- | ---- |
| [`processing-lab/`](processing-lab/) | C++ geometry, reflection prediction, integration, CONTROL execution and the CTest suite |
| [`mosflm/`](mosflm/) | Amalgamated Fortran, built only as the optional reference `ipmosflm` |
| `bin/ProcessingLab` | The processing executable (`ProcessingLab.exe` on Windows) |
| `bin/ipmosflm` | Optional Fortran reference (`make fortran-oracle`) |
| `processing-lab/build/processinglab_dev` | Developer diagnostics; not installed |

## Documentation

- [USAGE.md](USAGE.md): command-line options, run modes, metadata imports,
  CONTROL keywords and outputs.
- [DEVELOPMENT.md](DEVELOPMENT.md): direct CMake builds, the Fortran
  reference implementation and test fixtures.
- [processing-lab/README.md](processing-lab/README.md): the C++ library and
  `processinglab_dev`.

## License

BSD-style; see [COPYING](COPYING) for the terms and attribution
(A.G.W. Leslie / MRC LMB).
