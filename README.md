# ProcessingLab

Progressive C++ port of the MOSFLM diffraction-processing core, plus the
legacy Autotools `ipmosflm` tree used as the scientific oracle.

| Path | Role |
| ---- | ---- |
| [`mosflm-cpp/`](mosflm-cpp/) | C++ geometry, PREDICT, integration, C ABI, CTest |
| [`mosflm/`](mosflm/) | Amalgamated Fortran (oracle / Phase C splice) |
| [`porting.md`](porting.md) | Migration roadmap and current status |
| `bin/ipmosflm` | Legacy keyword-driven binary (Autotools) |

## Quick start (C++ port)

```bash
cmake -S mosflm-cpp -B mosflm-cpp/build -DCMAKE_BUILD_TYPE=Release
cmake --build mosflm-cpp/build
ctest --test-dir mosflm-cpp/build --output-on-failure
./mosflm-cpp/build/mosflm_cpp_cli version
```

Top-level CMake (`-S .`) adds `mosflm-cpp` when `PROCESSINGLAB_BUILD_CPP=ON`
(default).

## Legacy Autotools oracle

Requires CCP4. See [`porting.md`](porting.md) for macOS notes.

```bash
./configure --with-ccp4-prefix=/Applications/ccp4-9
make -C mosflm
# Experimental C++ stills PREDICT splice:
#   cmake --build mosflm-cpp/build --target mosflm_predict_hook
#   make -C mosflm MOSFLM_CPP_PREDICT=1 ipmosflm
```

## License

BSD-style — see [`COPYING`](COPYING) (A.G.W. Leslie / MRC LMB).
