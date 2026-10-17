# mosflm-cpp

Progressive C++ reimplementation of MOSFLM’s computational core. See
[`../porting.md`](../porting.md) for strategy and legacy build notes.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
./build/mosflm_cpp_cli
```

Requires CMake ≥ 3.24, a C++20 compiler, and network access on first configure
(FetchContent downloads Eigen 3.4.0).

## Layout

- `include/mosflm/geometry` — UnitCell, ReciprocalLattice, OrientationMatrix, DetectorGeometry
- `include/mosflm/crystal` — Gemmi-backed CrystalMetadata / SpaceGroup / MtzAdapter
- `include/mosflm/image` — DiffractionImage API stub (CBFlib reader TBD)
- `include/mosflm/legacy` — CELL(6) adapters / future Fortran bridge
- `tests/unit`, `tests/differential`
- `testdata/` — analytical fixtures, CBF sample symlink, oracle capture instructions
- `tools/cbf_smoke.c` — CBFlib header/pixel smoke test
- `docs/` — `FIRST_ROUTINE.md`, `CBF_STATUS.md`, `XDS_BRIDGE.md`
