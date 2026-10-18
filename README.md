# ProcessingLab

Progressive C++ port of the MOSFLM diffraction-processing core.

**Phase D Done:** default build/install is the C++ core. The release binary is
**`IntegrationLab`** (`IntegrationLab.exe` on Windows) — a thin CONTROL →
PROCESS driver that does **not** link the Fortran amalgamation. Autotools
`mosflm/` remains an optional scientific oracle named **`ipmosflm`** (F77).

| Path | Role |
| ---- | ---- |
| [`mosflm-cpp/`](mosflm-cpp/) | C++ geometry, PREDICT, integration, CONTROL, CTest |
| [`mosflm/`](mosflm/) | Amalgamated Fortran (**opt-in oracle** → `bin/ipmosflm`) |
| [`porting.md`](porting.md) | Migration roadmap and current status |
| `bin/IntegrationLab` | C++ CONTROL binary (default; `.exe` on Windows) |
| `bin/mosflm_cpp_cli` | Full C++ CLI |
| `bin/ipmosflm` | F77 oracle only (`make fortran-oracle`) |

## Quick start (default)

```bash
make cpp          # builds cbf/lib/libcbf.a then mosflm-cpp
make test         # ctest in mosflm-cpp/build
./bin/IntegrationLab version
./bin/mosflm_cpp_cli process-insulin
./bin/IntegrationLab HKLOUT /tmp/out.mtz --image /path/to/insulin_1_001.img \
  mosflm-cpp/testdata/images/known_indexed/cpp_process.in
```

If configuring `mosflm-cpp` directly, build CBFlib first: `make -C cbf`.

Top-level CMake (`-S .`) adds `mosflm-cpp` when `PROCESSINGLAB_BUILD_CPP=ON`
(default).

## Optional Fortran oracle

Requires CCP4. See [`porting.md`](porting.md) for macOS notes.

```bash
make fortran-oracle   # → bin/ipmosflm (F77; does not overwrite IntegrationLab)
# or:
#   ./configure --with-ccp4-prefix=/Applications/ccp4-9
#   make -C mosflm MOSFLM_FORTRAN_ORACLE=1 ipmosflm
# Experimental C++ PREDICT splice into the oracle:
#   cmake -S mosflm-cpp -B mosflm-cpp/build -DMOSFLM_FORTRAN_ORACLE=ON
#   cmake --build mosflm-cpp/build --target mosflm_predict_hook
#   make -C mosflm MOSFLM_FORTRAN_ORACLE=1 MOSFLM_CPP_PREDICT=1 ipmosflm
```

## License

BSD-style — see [`COPYING`](COPYING) (A.G.W. Leslie / MRC LMB).
