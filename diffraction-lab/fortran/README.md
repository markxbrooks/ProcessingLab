# Fortran ↔ C++ Phase C bridge

`mosflm_c_bind.f90` mirrors [`../include/mosflm/legacy/MosflmCApi.h`](../include/mosflm/legacy/MosflmCApi.h)
via `iso_c_binding` (`predict_still_cpp` / `predict_osc_cpp` convenience wrappers).

For fixed-form amalgamation callers, use the C shim `predict_still_cpp_`
(**gfortran** trailing-underscore ABI) — see `test_predict_still_f77.f`.
Amalgamation-oriented wrapper: `predict_cpp_still` in `predict_cpp_bridge.f90`
(see [`../docs/PHASE_C_SPLICE.md`](../docs/PHASE_C_SPLICE.md)).

ABI contract (status codes, count-only, UB/SG): [`../docs/C_ABI.md`](../docs/C_ABI.md).

GENDATA pack/unpack (mm ↔ 10 µm `XG`/`YG`): `predict_gendata_pack.f90`.  
Amalgamation hooks: `predict_cpp_hook.f90`, `predict_cpp_from_common.f90`
(`libmosflm_predict_hook` for Autotools). See [`../docs/PHASE_C_SPLICE.md`](../docs/PHASE_C_SPLICE.md).

```bash
cmake -S . -B build -DMOSFLM_WITH_FORTRAN=ON
cmake --build build --target mosflm_predict_hook
ctest --test-dir build -R 'predict_cpp_|ub_from_psi|c_api_contract|fortran_c_api' --output-on-failure
```

Requires `gfortran` (Homebrew GCC on macOS). Disable with `-DMOSFLM_WITH_FORTRAN=OFF`.
