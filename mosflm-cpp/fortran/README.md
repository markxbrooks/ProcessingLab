# Fortran ↔ C++ Phase C bridge

`mosflm_c_bind.f90` mirrors [`../include/mosflm/legacy/MosflmCApi.h`](../include/mosflm/legacy/MosflmCApi.h)
via `iso_c_binding` (`predict_still_cpp` / `predict_osc_cpp` convenience wrappers).

For fixed-form amalgamation callers, use the C shim `predict_still_cpp_`
(**gfortran** trailing-underscore ABI) — see `test_predict_still_f77.f`.
Amalgamation-oriented wrapper: `predict_cpp_still` in `predict_cpp_bridge.f90`
(see [`../docs/PHASE_C_SPLICE.md`](../docs/PHASE_C_SPLICE.md)).

ABI contract (status codes, count-only, UB/SG): [`../docs/C_ABI.md`](../docs/C_ABI.md).

```bash
cmake -S . -B build -DMOSFLM_WITH_FORTRAN=ON
cmake --build build --target test_fortran_c_api test_predict_still_f77 \
  test_predict_cpp_bridge test_predict_cpp_psi
ctest --test-dir build -R 'fortran_c_api|predict_still_f77|predict_cpp_bridge|predict_cpp_psi|ub_from_psi|c_api_contract' --output-on-failure
```

Requires `gfortran` (Homebrew GCC on macOS). Disable with `-DMOSFLM_WITH_FORTRAN=OFF`.
