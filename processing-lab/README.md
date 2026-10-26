# processing-lab (ProcessingLab)

C++ computational core for **ProcessingLab** (MOSFLM port). See
[`../DEVELOPMENT.md`](../DEVELOPMENT.md) for build, oracle and fixture notes.

**Status (Phase D Done / D6):** geometry + Gemmi crystal + miniCBF + **SMV** I/O
+ DPS Bravais façade + stills/oscillation `PREDICT` + **MEAS** boxes +
**Rossmann standard-profile** + Stage 2 **`PostRefiner`** + CONTROL
**`ControlSession`** + **`ProcessingLab`** (default release binary; no
amalgamation; `.exe` on Windows) + LP + C ABI + MTZ + pure-C++ `ProcessPipeline`.
F77 oracle remains `ipmosflm` (`make fortran-oracle`). Oracle splice library
`processinglab_predict_hook` only with `-DMOSFLM_FORTRAN_ORACLE=ON`. CTests include
`control_session`, `insulin_control`, `integrationlab`, `ipmosflm_perf`.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
# Optional CCTBX Python POC (Phenix cctbx.python):
#   cmake -S . -B build -DMOSFLM_WITH_CCTBX=ON
# Optional amalgamation splice archive (Autotools oracle only):
#   cmake -S . -B build -DMOSFLM_FORTRAN_ORACLE=ON
cmake --build build
ctest --test-dir build --output-on-failure
# Install ProcessingLab (processinglab_dev stays in build/):
#   cmake --install build --prefix /path/to/prefix
```

Requires CMake ≥ 3.24, a C++20 compiler, ZLIB, vendored `../third_party/cbf/lib/libcbf.a`
(`make -C ../cbf` or top-level `make cpp`). Eigen is taken from the system
(`find_package(Eigen3)`, e.g. Homebrew `eigen`) when present; otherwise a
3.4.0 tarball is downloaded (no git clone). Gemmi still uses FetchContent.
DPS Bravais tests need
`../lib/dps_index.a` and `../lib/libdps.a`. Optional **HDF5** (Homebrew /
CCP4) enables `Hdf5ImageReader` for Dectris Eiger NeXus masters
(`*_master.h5` / `.nxs`) via vendored bitshuffle+LZ4 filters; keep
`data_*.h5` beside the master. Soft-skip fixture test without
`MOSFLM_HDF5_FIXTURE`.

## CLI

`processinglab_dev` is a developer diagnostics tool (fixtures, smoke
checks); it is not installed. Use `ProcessingLab` to process data.

```bash
./build/processinglab_dev version
./build/processinglab_dev predict-still /tmp/still.dump
./build/processinglab_dev predict-osc /tmp/osc.dump
./build/processinglab_dev vertical-slice /tmp/slice.mtz
./build/processinglab_dev cbf-slice testdata/cbf/sample_minicbf.cbf /tmp/cbf.mtz
./build/processinglab_dev profile-fit /tmp/profile.mtz
./build/processinglab_dev profile-fit-cbf testdata/cbf/sample_minicbf.cbf /tmp/pf.mtz
# Requires CCP4 insulin SMV (or MOSFLM_IMAGE / CCP4 env):
./build/processinglab_dev process-insulin   # → insulin_cpp_process.mtz
./build/processinglab_dev process-insulin /path/to/insulin_1_001.img /tmp/out.mtz
# Default release binary (D6/E; ProcessingLab.exe on Windows):
./build/ProcessingLab HKLOUT /tmp/out.mtz --image /path/to/insulin_1_001.img \
  testdata/images/known_indexed/cpp_process.in
# Find spots → .spt → AUTOINDEX → multi-frame PROCESS (HDF5 master):
./build/ProcessingLab HKLOUT /tmp/out.mtz --image /path/to/master.h5 <<'EOF'
NPROC 8
FINDSPOTS THRESHOLD 5 FIND 1 FILE spots.spt MULTIFRAME RANGE 1 50
GO
AUTOINDEX DPS THRESHOLD 20
GO
PROCESS 1 TO 10
GO
EXIT
EOF
# OpenMP: brew install libomp (Apple Clang). See docs/XDS_BRIDGE.md § NPROC.
# Wall-time vs F77 ipmosflm (CCP4 default /apps/ccp4/ccp4-9/bin/ipmosflm):
./tools/bench_ipmosflm_perf.sh
# or: ctest --test-dir build -R ipmosflm_perf --output-on-failure
```

## Layout

| Path | Role |
| ---- | ---- |
| `include/mosflm/geometry` | UnitCell, ReciprocalLattice, Orientation, Detector, Dstar, Rotmat, Clcalc, Matset |
| `include/mosflm/crystal` | Gemmi CrystalMetadata / SpaceGroup / MtzAdapter |
| `include/mosflm/image` | `CbfImageReader`, `SmvImageReader`, `Hdf5ImageReader` (Eiger `*_master.h5` / `.nxs`) |
| `include/mosflm/diffraction` | Predict / integrate / `ProcessPipeline`; see `docs/DOMAIN_API.md` |
| `include/mosflm/refinement` | Stage 2 `PostRefiner` (POSTREFL-class φ LS) |
| `include/mosflm/control` | CONTROL keyword subset → `ProcessPipeline` |
| `include/mosflm/indexing` | BravaisClassifier + SUGGLATT select (penalty table log), `SpotFinder` + `.spt` I/O, CONTROL FINDSPOTS/AUTOINDEX |
| `include/mosflm/legacy` | CELL(6) + hardened `MosflmCApi.h` (see `docs/C_ABI.md`) |
| `fortran/` | `iso_c_binding` mirror + `predict_still_cpp` (Phase C; needs gfortran) |
| `tests/` | unit + differential (oracle dumps under `testdata/oracle/`) |
| `tools/` | `predict_*_oracle.py`, `cbf_smoke.c`, `cctbx_cell_poc.py` |
| `docs/` | FIRST_ROUTINE, CBF_STATUS, CCTBX_EVAL, XDS_BRIDGE |

## Bravais API

The Cambridge **5-argument** `bravais()` prototype lives in a single unique
header: [`../include/dps/bravais_api.h`](../include/dps/bravais_api.h).
Local `bravais.h` files under `util/`, `index/`, and `src/dps/bravais/` are
wrappers that include it (avoids divergent same-basename APIs).
