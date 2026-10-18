# Porting ProcessingLab (MOSFLM) to C++: a practical roadmap

**Repository:** ProcessingLab (formerly `ipmosflm`). Default `bin/ProcessingLab`
(`ProcessingLab.exe` on Windows) is the C++ CONTROL driver (Phase D / E); the
Autotools amalgamation is an optional F77 oracle (`make fortran-oracle` →
`bin/ipmosflm`). The C++ port lives in `mosflm-cpp/`.

I would approach this as a progressive reimplementation of MOSFLM's computational core, with a C++ interface and a compatibility layer around the existing crystallographic algorithms—not as a line-by-line translation of the entire repository.

The key distinction is that ProcessingLab is the project, while the `mosflm` directory is the part to investigate first for the underlying diffraction-processing algorithms. Before choosing an architecture, establish which source files implement the scientific calculations, which manage the application and its input/output, and which depend on legacy Fortran conventions.

Your first milestone should be a C++ executable that reproduces a small, independently verifiable MOSFLM calculation. That gives you a reliable foundation before you tackle image integration, refinement, or a graphical interface.

## Current status (Oct 2026)

Active work lives under [`mosflm-cpp/`](mosflm-cpp/). Build and run CTest there (optional `MOSFLM_WITH_CCTBX`).

| Topic | Status |
| ----- | ------ |
| C++ vertical slice | Geometry, miniCBF I/O, stills + oscillation `PREDICT`, box summation, 2D Gaussian profile fit, MTZ write — covered by CTest and `mosflm_cpp_cli` |
| PREDICT oracles | Analytical dumps + real insulin capture under `mosflm-cpp/testdata/oracle/known_indexed/`. CTest `insulin_hkl_overlap`: ASU-aware recall **100%** (`mosflm_lab_to_cpp`); `insulin_xdet`: median ~**23 px** (`ROTMAT2` + refined TILT/TWIST/YSCALE); `insulin_process_mtz`: exact-HKL recall ≥ 0.85 |
| Default `bin/ProcessingLab` | **C++** CONTROL → PROCESS (D6/E); no amalgamation link (`.exe` on Windows). F77 oracle: `make fortran-oracle` → `bin/ipmosflm`; optional `MOSFLM_CPP_PREDICT=1` splice ~89% exact-HKL (`known_indexed_cpp_predict/`) |
| Phase C | **Done** — lab-frame map, HKL/XDET differentials, mosaicity + absences in predictor, YSCALE/TILT/TWIST fields, PROCESS splice + MTZ CTest |
| Phase D / D1 | **Done** — pure C++ `ProcessPipeline` (SMV → predict → sum → MTZ); CTest `insulin_cpp_process` ASU recall ≥ 0.85 (no amalgamation hook) |
| Phase D / D2 | **Done** — `MeasIntegrator` (SEPARATION boxes, overlap, bg sigma-clip); CTest `insulin_isi` soft I/SIGI vs oracle; no Fortran MEAS |
| Phase D / D3 | **Done** — `StandardProfileFitter` (learned Rossmann profile); CTest `insulin_profile`; INTEG3 row = parity (global; not NMASKS) |

Phase table: see [§3](#3-choose-a-migration-strategy-that-allows-mixed-languages). Package layout and CLI: [`mosflm-cpp/README.md`](mosflm-cpp/README.md).

## Repo findings (verified locally)

| Item | Finding |
| ---- | ------- |
| Package | Repo **ProcessingLab**; Autotools still `AC_INIT([ipmosflm], [7.3.0])`; binary banner **Mosflm 7.4.0** (same as CCP4-9) |
| Core sources | Amalgamated Fortran in `mosflm/mosflm_all_ip_inc*.f` (~250k lines) + `unoptimizable.f` + `control.f` |
| Entry point | `PROGRAM MOSFLM` in `mosflm/unoptimizable.f:5`; command loop `SUBROUTINE CONTROL` in `control.f:4` |
| Indexing | DPS under `src/dps/` + `index/`; **REFIX stubbed** (`torefix*.f`) — public tree has no `refix_all.for` |
| GUI | External iMosflm (Tcl/Tk) over sockets; **no Tcl in this tree** |
| License | BSD-style `COPYING` (A.G.W. Leslie / MRC LMB) |
| C++ port home | `mosflm-cpp/` (this roadmap); top-level CMake project name `ProcessingLab` |

### High-value routine map (no-X amalgamation preferred)

| Routine | File | Approx. line | C++ status (Oct 2026) |
| ------- | ---- | ------------ | --------------------- |
| `PREDICT` | `mosflm_all_ip_inc_noX.f` | 171805 | Stills + oscillation + known-indexed HKL/XDET differentials; optional Fortran splice via `MOSFLM_CPP_PREDICT` |
| `DSTAR` | `unoptimizable.f` | 17692 | Thin helpers (`Dstar`) |
| `ROTMAT` | `unoptimizable.f` | 37480 | `rotmat_from_missetting` |
| `CLCALC` | `mosflm_all_ip_inc_noX.f` | 23322 | `cell_from_metric_tensor` |
| `MATSET` | `mosflm_all_ip_inc_noX.f` | 117930 | Thin `matset` (no COMMON derivatives) |
| `INTEG` / `INTEG2` | `unoptimizable.f` | 22256 / 20799 | Box summation integrator |
| `INTEG3` | `mosflm_all_ip_inc_noX.f` | 97457 | **Parity (global profile)** — `StandardProfileFitter` (Rossmann); Gaussian `ProfileFitter` kept for smoke tests |
| `MEAS` | `mosflm_all_ip_inc_noX.f` | 119120 | `MeasIntegrator` (SEPARATION boxes / overlap / bg clip; not full stripe buffer) |
| `PROCESS` | `unoptimizable.f` | 25854 | `ProcessPipeline` (predict → profile → MTZ; optional `PostRefiner`) |
| `POSTREFL` | amalgamation | ~168264 | **Stage 2** — `PostRefiner` (misset + cell scale + mosaic from φ; off by default in PROCESS) |
| `CONTROL` | `control.f` / amalgamation | — | **D5 subset** — `ControlSession` (CELL…EXIT → PROCESS); CLI `run` |
| `METRIC` / `INVCEL` | `torefix_noX.f` | 182 / 130 | `UnitCell` (+ differential test) |

Modular C oracles (good early differential targets): `src/dps/util/vector_math.c`, `cryst_matrix.c`, `mat_inversion.c`.

### Bravais header hygiene

There are multiple `bravais.h` basenames (`index/`, `util/`, `src/dps/bravais/`). The **definition** in `index/bravais.c` is the **5-argument** Cambridge API. Those local headers are now **wrappers** around a single prototype in [`include/dps/bravais_api.h`](include/dps/bravais_api.h). `bravais_fortran.c` calls the 5-arg API (with an explicit cast for Fortran’s flat `sys_id` buffer) and sets `*error = 0`.

### Building the legacy oracle (macOS arm64, Oct 2026)

```bash
export CCP4=/Applications/ccp4-9
export CFLAGS="-g -O0 -Wno-implicit-function-declaration -Wno-return-type"
export DYLD_LIBRARY_PATH="$CCP4/lib:/opt/homebrew/opt/hdf5/lib"

cd /path/to/ipmosflm
./configure --with-ccp4-prefix="$CCP4"

# Bundled JPEG (generates jpg/jconfig.h)
( cd jpg && ./configure --disable-shared && make libjpeg.a )

# Homebrew HDF5 headers/libs required for fromhdf5.c
# mosflm/Makefile was adjusted for: -I../jpg, -I/opt/homebrew/opt/hdf5/include,
# -lhdf5_hl -lhdf5, ../jpg/libjpeg.a, and to drop broken -lto_library.

make -j2
# → bin/ipmosflm

DYLD_LIBRARY_PATH="$CCP4/lib:/opt/homebrew/opt/hdf5/lib" ./bin/ipmosflm <<'EOF'
EXIT
EOF
```

Smoke-tested: version banner prints and `EXIT` terminates cleanly. Keyword smoke (`CELL` / `MATRIX` / Pilatus template errors) matches CCP4-9’s `ipmosflm`.

### Building mosflm-cpp (ProcessingLab)

```bash
# From the ProcessingLab repo root (or mosflm-cpp/):
cmake -S mosflm-cpp -B mosflm-cpp/build -DCMAKE_BUILD_TYPE=Release
# Optional: -DMOSFLM_WITH_CCTBX=ON  (Phenix cctbx.python POC)
cmake --build mosflm-cpp/build
ctest --test-dir mosflm-cpp/build --output-on-failure
```

After a directory rename, delete `mosflm-cpp/build` once so CMake does not keep
stale absolute paths from the old tree.

## 1. Start with a source-code inventory

Repository: [ipmosflm — mosflm source directory](https://ccp4forge.rc-harwell.ac.uk/ccp4/ipmosflm/-/tree/main/mosflm).

Don't begin by translating files in alphabetical order. First establish the language mix, build system, dependency graph, and execution flow.

Build a table with these categories:

| Category                  | What to identify                                            | Porting priority |
| ------------------------- | ----------------------------------------------------------- | ---------------- |
| Numerical kernels         | Geometry, coordinate transforms, least-squares calculations | High             |
| Diffraction processing    | Spot prediction, profile fitting, integration               | High             |
| Refinement                | Beam, detector, crystal, and orientation parameters         | High             |
| Data and I/O              | Reflection tables, images, experiment metadata              | Medium           |
| Runtime and orchestration | Command interpreter, global state, legacy control flow      | Medium           |
| User interface            | Tcl/Tk, scripts, display logic                              | Low initially    |

## 2. Preserve a working scientific reference

Before changing the implementation, get the original program to build and run. Keep a known-good executable and a small set of representative diffraction datasets.

MOSFLM's documented workflow has three major stages:

**Stage 1 — Autoindexing** — lattice, unit cell, orientation from spots.

**Stage 2 — Post-refinement** — **Done** — `PostRefiner` (cell scale, missetting, mosaic from φ residuals); CTests `post_refinement`, `insulin_postref` (not default-on in PROCESS).

**Stage 3 — Integration** — predict, measure intensities, write reflections.

For a C++ port, use **integration of a known, already indexed dataset** as the first substantial end-to-end target. Autoindexing adds ambiguity that obscures early discrepancies.

Regression corpus lives under [`mosflm-cpp/testdata/`](mosflm-cpp/testdata/README.md). Analytical oracles are in place; a known-indexed image set is still required for Fortran differentials (`tools/capture_ipmosflm_oracle.sh`).

## 3. Choose a migration strategy that allows mixed languages

Staged migration with **current status**:

| Phase | Intent | Status (Oct 2026) |
| ----- | ------ | ----------------- |
| **A** | Hybrid CMake build; narrow Fortran/C++ adapters | **Done** — `mosflm_geometry`, `LegacyMosflmAdapter` (CELL(6), predict, summation, profile), optional `mosflm_indexing` |
| **B** | Replace isolated kernels; keep Fortran as oracle | **Done** (advanced) — geometry, Gemmi crystal, miniCBF, stills+osc `PREDICT`, summation + profile fit, MTZ, Bravais façade, Stage 2 `PostRefiner`. Domain API: `docs/DOMAIN_API.md` |
| **C** | Switch implementation behind interfaces after differential tests pass | **Done** — `MOSFLM_CPP_PREDICT` splices GENERATE + REFLFIND; insulin PROCESS ~89% exact-HKL recall; CTests `insulin_hkl_overlap` / `insulin_xdet` / `insulin_process_mtz`. REEK divergence + ``ROTMAT2`` are in Phase D (`docs/PHASE_C_SPLICE.md`) |
| **D** | Retire Fortran only when all callers migrated | **Done** — D0–D6; default release is C++ core + thin `ipmosflm`; amalgamation opt-in oracle |

Keep COMMON blocks, column-major layout, and mangling inside `include/mosflm/legacy/`.

## 4. Design the C++ domain model before porting algorithms

Layout under `mosflm-cpp/` (populated): `geometry`, `diffraction`, `crystal`, `image`, `indexing`, `refinement`, `control`, `legacy`; CMake + CTest; Eigen + Gemmi via FetchContent. Optional CCTBX remains a standalone evaluator (`docs/CCTBX_EVAL.md`), not a hard dependency.

## 5. Build a vertical slice, not a collection of translated functions

| Step | Target | Status |
| ---- | ------ | ------ |
| 1 | Read diffraction image + metadata | **Done** — `CbfImageReader` (CBFlib + miniCBF `#` comments) |
| 2 | Experimental geometry | **Done** — `UnitCell` / `OrientationMatrix` / `DetectorGeometry` (+ Dstar, Rotmat, Clcalc, Matset) |
| 3 | Predict reflections (`PREDICT`) | **Done** — stills + oscillation; analytical oracles; see [`FIRST_ROUTINE.md`](mosflm-cpp/docs/FIRST_ROUTINE.md) |
| 4 | Summation integrate | **Done** — `SummationIntegrator` + `MeasIntegrator`. Profile: **Done** — `StandardProfileFitter` (INTEG3 Rossmann global profile) |
| 5 | Write MTZ; compare with numeric tolerances | **Done** — synthetic / CBF smoke + CTest `insulin_process_mtz` (cpp_predict vs Fortran oracle, exact-HKL recall ≥ 0.85) |

## First Fortran routine to replace: `PREDICT`

See [mosflm-cpp/docs/FIRST_ROUTINE.md](mosflm-cpp/docs/FIRST_ROUTINE.md) for callers, I/O contract, and regression criteria. C++ stills/oscillation path, dump format, and Phase C Fortran splice (`MOSFLM_CPP_PREDICT`) are implemented.

## Milestones

### First week (done)

Days 1–2: Source audit + baseline build.

Days 2–3: Regression corpus scaffolding.

Days 3–4: First C++ geometry kernel + unit tests.

Day 5: Document first legacy routine to replace (`PREDICT`).

### Next

1. Real known-indexed images — **done** (`capture_ipmosflm_oracle.sh` → `testdata/oracle/known_indexed/`).
2. **Phase C** — **Done** (D0 closed): HKL/XDET/MTZ CTests, mosaicity + absences, YSCALE pixel aspect + TILT/TWIST fields, PROCESS splice ~89% exact-HKL recall.
3. **Phase D / D1–D6** — **Done**: SMV + `ProcessPipeline` + MEAS + Rossmann profile + Stage 2 `PostRefiner` + CONTROL + `ProcessingLab` / `make cpp`; amalgamation behind `MOSFLM_FORTRAN_ORACLE` → `ipmosflm`. CTests include `insulin_control`, `integrationlab`. See [`portingD.md`](portingD.md).
4. **Phase E** — CCP4 drop-in roadmap: see [`portingE.md`](portingE.md) (E0 CLI UX through E7 packaging).
