# mosflm-cpp (ProcessingLab)

C++ computational core for **ProcessingLab** (MOSFLM port). See
[`../porting.md`](../porting.md) for strategy and legacy build notes.

**Status (Phase D / D1 Done):** geometry + Gemmi crystal + miniCBF + **SMV** I/O
+ DPS Bravais façade + stills/oscillation `PREDICT` (mosaicity / absences) +
box summation + **2D Gaussian profile fit** + LP correction + C ABI + MTZ write
+ amalgamation `MOSFLM_CPP_PREDICT` hook + pure-C++ `ProcessPipeline`
(`process-insulin`). Known-indexed CTests: `insulin_hkl_overlap`,
`insulin_xdet`, `insulin_process_mtz`, `insulin_cpp_process`. Next: D2 MEAS
(`../portingD.md`).

Phase A–D table and vertical-slice checklist: [`../porting.md`](../porting.md)
(section *Current status* and §3 / §5).

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
# Optional CCTBX Python POC (Phenix cctbx.python):
#   cmake -S . -B build -DMOSFLM_WITH_CCTBX=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

Requires CMake ≥ 3.24, a C++20 compiler, ZLIB, vendored `../cbf/lib/libcbf.a`,
and network access on first configure (FetchContent: Eigen, Gemmi). DPS Bravais
tests need `../lib/dps_index.a` and `../lib/libdps.a`.

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
```

## Layout

| Path | Role |
| ---- | ---- |
| `include/mosflm/geometry` | UnitCell, ReciprocalLattice, Orientation, Detector, Dstar, Rotmat, Clcalc, Matset |
| `include/mosflm/crystal` | Gemmi CrystalMetadata / SpaceGroup / MtzAdapter |
| `include/mosflm/image` | `CbfImageReader`, `SmvImageReader` (ADSC insulin) |
| `include/mosflm/diffraction` | Predict / integrate / `ProcessPipeline`; see `docs/DOMAIN_API.md` |
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
