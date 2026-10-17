# Phase C: splicing C++ stills PREDICT into the amalgamation

## Goal

Route one stills-oriented call site from Fortran `PREDICT` through the C++ ABI
without rewriting the full reflection-list machinery in COMMON.

## Hook location

In `mosflm/mosflm_all_ip_inc_noX.f` near line **194120**:

```fortran
        CALL PREDICT(PSI, IFLAGPR)
```

## Blockers before applying the ifdef

The call site at ≈194120 sits inside interactive re-predict (fills COMMON
`/IOSP/` `NSPOT` and spot lists used by XML / pixel matching below). The C++
bridge returns flat `H,K,L,XMM,…` arrays; it does **not** yet:

1. Pack detector / cell / φ / `DMIN` from COMMON into bridge scalars
2. Scatter bridge outputs into the MOSFLM reflection COMMON arrays
3. Link `libmosflm_geometry` + `predict_cpp_bridge.f90` into Autotools `ipmosflm`

Until (1)–(3) exist, leave the amalgamation on `CALL PREDICT`. The ifdef
recipe below is ready to paste when packing is implemented.

### COMMON pack / unpack sketch

Fortran `PREDICT` does not compute spots itself; it calls `REEK` (then
`OVERLAP`), which mutates shared state. A C++ splice must mimic the **outputs
consumers read**, not the REEK call tree.

| Symbol | Role after `PREDICT` | C++ bridge source |
| ------ | -------------------- | ----------------- |
| `NSPOT` (`COMMON /IOSP/`) | Spot count | `nspot` |
| `IFLAGPR` | −1 if list overflow | `ierr != 0` or buffer truncate |
| Spot miller / detector lists (`XG`/`YG`, indices, φ beg/end, full flags — names vary by include) | Per-spot geometry | `h,k,l,x_mm,y_mm,phi_beg,phi_end,ifull` |

**Inputs to pack from COMMON** (typical): `CELL(6)`, detector distance /
wavelength / pixel / beam, still `PHI`, resolution / Ewald limits, and `PSI(3)`
(already the `PREDICT` argument). Optional: precomputed `AMAT`/`MATRIX` via
`mosflm_ub_from_matrix9` instead of PSI.

**Suggested implementation order:** (a) Fortran helper that copies COMMON →
bridge args and bridge → `NSPOT` + spot arrays; (b) unit-test that helper
against analytical stills fixture offline; (c) then enable
`MOSFLM_CPP_PREDICT` at ≈194120.

## Proposed ifdef (not applied to the Autotools tree yet)

```fortran
#ifdef MOSFLM_CPP_PREDICT
C---- Experimental stills PREDICT via mosflm-cpp (PSI → UB).
C     Pack CELL / detector from COMMON into locals, then:
        CALL PREDICT_CPP_STILL_PSI(CELL, DIST, WAVE, PX, PY, NF, NS,
     +       BX, BY, PHI, DMIN, ETOL, PSI, H, K, L, XMM, YMM, PHIB,
     +       PHIE, IFULL, NMAX, NSPOT, IERR)
        IF (IERR .NE. 0) IFLAGPR = -1
#else
        CALL PREDICT(PSI, IFLAGPR)
#endif
```

Bridge source: [`../fortran/predict_cpp_bridge.f90`](../fortran/predict_cpp_bridge.f90)
(`predict_cpp_still` identity UB; `predict_cpp_still_psi` with missetting).  
C ABI: `predict_still_cpp_` / `predict_still_psi_cpp_` /
`mosflm_ub_from_psi` / `mosflm_ub_from_matrix9` in
[`../include/mosflm/legacy/MosflmCApi.h`](../include/mosflm/legacy/MosflmCApi.h).

## Mapping notes (COMMON → bridge args)

| Bridge arg | Typical MOSFLM source | Notes |
| ---------- | --------------------- | ----- |
| `CELL(6)` | Direct cell COMMON | `REAL` Å / degrees |
| `DIST`, `WAVE`, `PX`, `PY` | Detector / beam keywords | mm, Å |
| `NF`, `NS` | Image dimensions | fast / slow |
| `BX`, `BY` | Beam centre | **beam-relative mm** in C++ (often 0,0 in that frame) |
| `PHI` | Still φ / image start | degrees |
| `PSI(3)` | Missetting at `PREDICT` | degrees → `UB = R(PSI)·A*` via `mosflm_ub_from_psi` |
| `DMIN`, `ETOL` | Resolution / Ewald | Å, Å⁻¹ |
| `H,K,L,XMM,…` | Output lists | C++ writes up to `NMAX` |

**MATRIX path:** if COMMON already holds a lab A/UB matrix (9 reals), copy with
`mosflm_ub_from_matrix9` and call `mosflm_predict_still` directly (no PSI).
Zero PSI matches the identity-UB bridge (`predict_cpp_still`).

## Linkage

The Autotools `ipmosflm` binary does **not** yet link `libmosflm_geometry.a`.
Steps when enabling the hook:

1. Build `mosflm-cpp` → `libmosflm_geometry.a`
2. Compile `predict_cpp_bridge.f90` with the amalgamation objects
3. Define `-DMOSFLM_CPP_PREDICT` and link C++ stdlib + Eigen/Gemmi deps as needed
4. Keep the `#else` path as the production default until differentials pass

## Verification in mosflm-cpp

```bash
ctest --test-dir build -R 'predict_cpp_bridge|predict_cpp_psi|ub_from_psi|predict_still_f77|fortran_c_api'
```
