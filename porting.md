# Porting ipmosflm to C++: a practical roadmap

I would approach this as a progressive reimplementation of MOSFLM's computational core, with a C++ interface and a compatibility layer around the existing crystallographic algorithms—not as a line-by-line translation of the entire repository.

The key distinction is that `ipmosflm` is the project, while the `mosflm` directory is the part to investigate first for the underlying diffraction-processing algorithms. Before choosing an architecture, establish which source files implement the scientific calculations, which manage the application and its input/output, and which depend on legacy Fortran conventions.

Your first milestone should be a C++ executable that reproduces a small, independently verifiable MOSFLM calculation. That gives you a reliable foundation before you tackle image integration, refinement, or a graphical interface.

## Repo findings (verified locally)

| Item | Finding |
| ---- | ------- |
| Package | Autotools **ipmosflm 7.3.0** (`configure.ac`); binary `bin/ipmosflm` |
| Core sources | Amalgamated Fortran in `mosflm/mosflm_all_ip_inc*.f` (~250k lines) + `unoptimizable.f` + `control.f` |
| Entry point | `PROGRAM MOSFLM` in `mosflm/unoptimizable.f:5`; command loop `SUBROUTINE CONTROL` in `control.f:4` |
| Indexing | DPS under `src/dps/` + `index/`; **REFIX stubbed** (`torefix*.f`) — public tree has no `refix_all.for` |
| GUI | External iMosflm (Tcl/Tk) over sockets; **no Tcl in this tree** |
| License | BSD-style `COPYING` (A.G.W. Leslie / MRC LMB) |
| C++ port home | `mosflm-cpp/` (this roadmap) |

### High-value routine map (no-X amalgamation preferred)

| Routine | File | Approx. line |
| ------- | ---- | ------------ |
| `PREDICT` | `mosflm_all_ip_inc_noX.f` | 171805 |
| `DSTAR` | `unoptimizable.f` | 17692 |
| `ROTMAT` | `unoptimizable.f` | 37480 |
| `CLCALC` | `mosflm_all_ip_inc_noX.f` | 23322 |
| `MATSET` | `mosflm_all_ip_inc_noX.f` | 117930 |
| `INTEG` / `INTEG2` | `unoptimizable.f` | 22256 / 20799 |
| `INTEG3` | `mosflm_all_ip_inc_noX.f` | 97457 |
| `MEAS` | `mosflm_all_ip_inc_noX.f` | 119120 |
| `PROCESS` | `unoptimizable.f` | 25854 |
| `METRIC` / `INVCEL` | `torefix_noX.f` | 182 / 130 |

Modular C oracles (good early differential targets): `src/dps/util/vector_math.c`, `cryst_matrix.c`, `mat_inversion.c`.

### Bravais header hygiene

There are multiple `bravais.h` basenames (`index/`, `util/`, `src/dps/bravais/`). The **definition** in `index/bravais.c` is the **5-argument** Cambridge API. `util/bravais.h` previously declared a stale 6-arg form; `bravais_fortran.c` was updated to call the 5-arg API and set `*error = 0` for Fortran. Prefer a single include path under `src/dps/bravais/` long-term rather than divergent same-basename headers.

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

Smoke-tested: version banner prints and `EXIT` terminates cleanly.

### Building mosflm-cpp

```bash
cd mosflm-cpp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

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

**Stage 2 — Post-refinement** — cell, orientation, mosaic spread from partials.

**Stage 3 — Integration** — predict, measure intensities, write reflections.

For a C++ port, use **integration of a known, already indexed dataset** as the first substantial end-to-end target. Autoindexing adds ambiguity that obscures early discrepancies.

Regression corpus lives under [`mosflm-cpp/testdata/`](mosflm-cpp/testdata/README.md).

## 3. Choose a migration strategy that allows mixed languages

Staged migration:

- **Phase A** — Hybrid CMake build; narrow Fortran/C++ adapters.
- **Phase B** — Replace isolated kernels; keep Fortran as oracle.
- **Phase C** — Switch implementation behind interfaces after differential tests pass.
- **Phase D** — Retire Fortran only when all callers migrated.

Keep COMMON blocks, column-major layout, and mangling inside `include/mosflm/legacy/`.

## 4. Design the C++ domain model before porting algorithms

Layout under `mosflm-cpp/` (created): geometry, diffraction, refinement, image, io, legacy adapters; CMake + CTest; Eigen via FetchContent.

## 5. Build a vertical slice, not a collection of translated functions

1. Read diffraction image + metadata.
2. Represent experimental geometry (this milestone: UnitCell / Orientation / Detector).
3. Predict reflections (`PREDICT` path — see below).
4. Summation integrate (defer profile fitting).
5. Write MTZ; compare with numeric tolerances.

## First Fortran routine to replace: `PREDICT`

See [mosflm-cpp/docs/FIRST_ROUTINE.md](mosflm-cpp/docs/FIRST_ROUTINE.md) for callers, I/O contract, and regression criteria.

## My recommendation for your first week

Days 1–2: Source audit + baseline build (done).

Days 2–3: Regression corpus scaffolding (done).

Days 3–4: First C++ geometry kernel + unit tests (done).

Day 5: Document first legacy routine to replace (`PREDICT`) (done).
