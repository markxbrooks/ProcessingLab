# Developing ProcessingLab

Notes for working on ProcessingLab itself. For building and running it,
see [README.md](README.md). For the command-line reference, see
[USAGE.md](USAGE.md).

## Direct CMake builds

The top-level `make` drives CMake in `processing-lab/build`. To configure
`processing-lab` yourself, build the vendored CBFlib first:

```bash
make -C cbf
cmake -S processing-lab -B processing-lab/build -DCMAKE_BUILD_TYPE=Release
cmake --build processing-lab/build
ctest --test-dir processing-lab/build --output-on-failure
```

- **Top-level CMake:** `-S .` adds `processing-lab` when
  `PROCESSINGLAB_BUILD_CPP=ON` (the default).
- **Dependencies:** Eigen comes from the system when present; otherwise a
  3.4.0 tarball is downloaded. Gemmi uses FetchContent.
- **DPS:** the Bravais and AutoIndex tests need `lib/dps_index.a` and
  `lib/libdps.a`, from `make dps-lib` (gfortran).
- **HDF5:** HDF5 from Homebrew or CCP4 enables Eiger/NeXus masters.
- **Test output:** tests write under `processing-lab/build/test_output/`,
  never into the source tree.
- **`processinglab_dev`:** developer diagnostics (analytical prediction
  dumps, smoke slices), built in `processing-lab/build/` and never
  installed. See [processing-lab/README.md](processing-lab/README.md).

## Fortran reference implementation (`ipmosflm`)

The amalgamated Fortran under `mosflm/` is built only as a reference for
differential tests. It needs CCP4, and doesn't replace `ProcessingLab`:

```bash
make fortran-oracle                  # → bin/ipmosflm
# or, configured against a CCP4 tree:
./configure --with-ccp4-prefix=/Applications/ccp4-9
make -C mosflm MOSFLM_FORTRAN_ORACLE=1 ipmosflm
```

To compare the wall time of `ProcessingLab` and `ipmosflm` on a
known-indexed PROCESS (CCP4's `ipmosflm` by default):

```bash
./processing-lab/tools/bench_ipmosflm_perf.sh
```

An experimental build splices the C++ `PREDICT` into the Fortran program:

```bash
cmake -S processing-lab -B processing-lab/build -DMOSFLM_FORTRAN_ORACLE=ON
cmake --build processing-lab/build --target processinglab_predict_hook
make -C mosflm MOSFLM_FORTRAN_ORACLE=1 MOSFLM_CPP_PREDICT=1 ipmosflm
```

## Test fixtures

Some differential tests compare against `ipmosflm` runs on the CCP4 insulin
example. Their MTZs aren't in git. Regenerate them with
`processing-lab/tools/capture_ipmosflm_oracle.sh`. On macOS the script sets
the library path itself, and if a step fails it prints the end of the
`ipmosflm` log.

- **`testdata/oracle/known_indexed/insulin_1_001.mtz`:** run the script
  with CCP4's `ipmosflm`, the default when `bin/ipmosflm` doesn't exist.
- **`testdata/oracle/known_indexed_cpp_predict/insulin_1_001.mtz`:** build
  `ipmosflm` with `MOSFLM_CPP_PREDICT=1` (above), run the script with
  `IPMOSFLM=` pointing at it, and copy the MTZ there.

Tests that need data which isn't present (the insulin image, the thaumatin
master, the 76c1b sweep) skip themselves.
