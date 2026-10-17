# Fortran ↔ C++ Phase C bridge

`mosflm_c_bind.f90` mirrors [`../include/mosflm/legacy/MosflmCApi.h`](../include/mosflm/legacy/MosflmCApi.h)
via `iso_c_binding`. `predict_still_cpp` is the drop-in style stills PREDICT
entry that future amalgamation call sites can route through.

```bash
cmake -S . -B build -DMOSFLM_WITH_FORTRAN=ON
cmake --build build --target test_fortran_c_api
ctest --test-dir build -R fortran_c_api --output-on-failure
```

Requires `gfortran` (Homebrew GCC on macOS). Disable with `-DMOSFLM_WITH_FORTRAN=OFF`.
