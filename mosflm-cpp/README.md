# mosflm-cpp (ProcessingLab)

C++ computational core for **ProcessingLab** (MOSFLM port). See
[`../porting.md`](../porting.md) for strategy and legacy build notes.

**Status (Phase D Done / D6):** geometry + Gemmi crystal + miniCBF + **SMV** I/O
+ DPS Bravais façade + stills/oscillation `PREDICT` + **MEAS** boxes +
**Rossmann standard-profile** + Stage 2 **`PostRefiner`** + CONTROL
**`ControlSession`** + **`IntegrationLab`** (default release binary; no
amalgamation; `.exe` on Windows) + LP + C ABI + MTZ + pure-C++ `ProcessPipeline`.
F77 oracle remains `ipmosflm` (`make fortran-oracle`). Oracle splice library
`mosflm_predict_hook` only with `-DMOSFLM_FORTRAN_ORACLE=ON`. CTests include
`control_session`, `insulin_control`, `integrationlab`.

Phase A–D table and vertical-slice checklist: [`../porting.md`](../porting.md)
(section *Current status* and §3 / §5).

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
# Optional CCTBX Python POC (Phenix cctbx.python):
#   cmake -S . -B build -DMOSFLM_WITH_CCTBX=ON
# Optional amalgamation splice archive (Autotools oracle only):
#   cmake -S . -B build -DMOSFLM_FORTRAN_ORACLE=ON
cmake --build build
ctest --test-dir build --output-on-failure
# Install IntegrationLab + mosflm_cpp_cli:
#   cmake --install build --prefix /path/to/prefix
```

Requires CMake ≥ 3.24, a C++20 compiler, ZLIB, vendored `../cbf/lib/libcbf.a`
(`make -C ../cbf` or top-level `make cpp`), and network access on first
configure (FetchContent: Eigen, Gemmi). DPS Bravais tests need
`../lib/dps_index.a` and `../lib/libdps.a`.

## CLI

```bash
./build/mosflm_cpp_cli version
./build/mosflm_cpp_cli predict-still /tmp/still.dump
./build/mosflm_cpp_cli predict-osc /tmp/osc.dump
./build/mosflm_cpp_cli vertical-slice /tmp/slice.mtz
./build/mosflm_cpp_cli cbf-slice testdata/cbf/sample_minicbf.cbf /tmp/cbf.mtz
./build/mosflm_cpp_cli profile-fit /tmp/profile.mtz
./build/mosflm_cpp_cli profile-fit-cbf testdata/cbf/sample_minicbf.cbf /tmp/pf.mtz
# Requires CCP4 insulin SMV (or MOSFLM_IMAGE / CCP4 env):
./build/mosflm_cpp_cli process-insulin   # → insulin_cpp_process.mtz
./build/mosflm_cpp_cli process-insulin /path/to/insulin_1_001.img /tmp/out.mtz
# Keyword-driven PROCESS (D5):
./build/mosflm_cpp_cli run testdata/images/known_indexed/cpp_process.in \
  --image /path/to/insulin_1_001.img /tmp/control.mtz
# Default release binary (D6/E; IntegrationLab.exe on Windows):
./build/IntegrationLab HKLOUT /tmp/out.mtz --image /path/to/insulin_1_001.img \
  testdata/images/known_indexed/cpp_process.in
```

## Layout

| Path | Role |
| ---- | ---- |
| `include/mosflm/geometry` | UnitCell, ReciprocalLattice, Orientation, Detector, Dstar, Rotmat, Clcalc, Matset |
| `include/mosflm/crystal` | Gemmi CrystalMetadata / SpaceGroup / MtzAdapter |
| `include/mosflm/image` | `CbfImageReader`, `SmvImageReader` (ADSC insulin) |
| `include/mosflm/diffraction` | Predict / integrate / `ProcessPipeline`; see `docs/DOMAIN_API.md` |
| `include/mosflm/refinement` | Stage 2 `PostRefiner` (POSTREFL-class φ LS) |
| `include/mosflm/control` | CONTROL keyword subset → `ProcessPipeline` |
| `include/mosflm/indexing` | BravaisClassifier (DPS façade) |
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
