# Mosflm C ABI contract

Public header: [`../include/mosflm/legacy/MosflmCApi.h`](../include/mosflm/legacy/MosflmCApi.h).

The C++ domain layer uses typed request objects
([`DOMAIN_API.md`](DOMAIN_API.md)); this ABI stays flat and adapts into those
requests.

## Status codes

| Code | Name | Meaning |
| ---- | ---- | ------- |
| 0 | `MOSFLM_SUCCESS` | Completed; see `out_n` |
| 1 | `MOSFLM_ERROR_INVALID_ARGUMENT` | Null/invalid inputs or non-finite geometry |
| 2 | `MOSFLM_ERROR_BUFFER_TOO_SMALL` | Truncated write; `out_n` = elements written |
| 3 | `MOSFLM_ERROR_INTERNAL` | Caught C++ exception at the ABI boundary |

## Output counts

- `out_n` is always the **number written** into the caller buffer (0 on error when writable).
- It is zeroed as soon as `out_n` is non-null, before further validation.
- **Count-only** predict: `out_cap == 0` and `out == NULL` → `out_n` = total spots, return `SUCCESS`.

## Orientation and space group

- Pass `mosflm_ub9*` (row-major lab UB, `r* = UB * h`) or `NULL` for the canonical cell orientation (`R = I`).
- Build UB from MOSFLM missetting: `mosflm_ub_from_psi(cell, psi_deg[3], &ub)` → `UB = R(PSI)·A*` (`ROTMAT`).
- Copy an existing MATRIX: `mosflm_ub_from_matrix9(m9, &ub)`.
- F77 with PSI: `predict_still_psi_cpp_` (zero PSI ≡ `predict_still_cpp_`).
- Spindle φ is applied about +X on top of UB.
- `spacegroup` `NULL`/empty → `"P 1"`. Systematic absences are **not** applied yet.

## LP correction

Domain `LpCorrection` is the single source of truth. `mosflm_lp_factor` delegates to it. Near-beam `|sin 2θ|` clamp returns `1.0` (not physically valid on-axis).

## F77 shim

`predict_still_cpp_` (gfortran trailing underscore) is production-only: caller supplies `CELL(6)` and detector scalars. Fill inputs via `mosflm_analytical_still_fixture` from C/tests, or hard-code geometry as in `fortran/test_predict_still_f77.f`.
