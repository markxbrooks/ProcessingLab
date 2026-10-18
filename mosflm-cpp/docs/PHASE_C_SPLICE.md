# Phase C: splicing C++ stills PREDICT into the amalgamation

## Goal

Route one stills-oriented call site from Fortran `PREDICT` through the C++ ABI
without rewriting the full reflection-list machinery in COMMON.

## Hook location (applied)

In `mosflm/mosflm_all_ip_inc_noX.f` near line **194120** (`REFLFIND`), behind
`#ifdef MOSFLM_CPP_PREDICT` (default **off** — production still calls `PREDICT`).

`REFLFIND` does not include `/CELLCOM/` or `/XY/`. The ifdef calls
`PREDICT_CPP_FROM_COMMON`, which reads `CELL` / `XTOFD` from those commons and
writes into the in-scope `/GENDATA/` + `/IOSP/` arrays.

```fortran
#ifdef MOSFLM_CPP_PREDICT
        DMINCPP = 0.0
        IF (DSTMAX.GT.1.0E-12) DMINCPP = WAVE/DSTMAX
        ETOLCPP = 0.0025
        CALL PREDICT_CPP_FROM_COMMON(PSI, IFLAGPR, WAVE, RAST,
     +       RASTY, IYLEN, NREC, PHIBEG, DMINCPP, ETOLCPP, NREFLS,
     +       IHG, IKG, ILG, XG, YG, PHIG, PHIWG, IMG, NSPOT)
#else
        CALL PREDICT(PSI, IFLAGPR)
#endif
```

## Stack

| Layer | File |
| ----- | ---- |
| Amalgamation ifdef | `mosflm/mosflm_all_ip_inc_noX.f` |
| COMMON reader | [`../fortran/predict_cpp_from_common.f90`](../fortran/predict_cpp_from_common.f90) |
| Hook | [`../fortran/predict_cpp_hook.f90`](../fortran/predict_cpp_hook.f90) |
| GENDATA pack | [`../fortran/predict_gendata_pack.f90`](../fortran/predict_gendata_pack.f90) |
| C ABI | [`../include/mosflm/legacy/MosflmCApi.h`](../include/mosflm/legacy/MosflmCApi.h) |

## Enabling the experimental path

Verified on macOS arm64 (Oct 2026): binary links and starts; `nm` shows
`_predict_cpp_from_common_` / `_predict_still_psi_cpp_`.

```bash
# 1. Build C++ + Fortran hook archive
cmake -S mosflm-cpp -B mosflm-cpp/build -DMOSFLM_WITH_FORTRAN=ON
cmake --build mosflm-cpp/build --target mosflm_predict_hook

# 2. Rebuild amalgamation with -DMOSFLM_CPP_PREDICT and link hooks
cd mosflm
rm -f mosflm_all_ip_inc_noX.o
make MOSFLM_CPP_PREDICT=1 ipmosflm
```

`mosflm/Makefile` always passes `-cpp` so `#ifdef` is safe when the macro is
undefined. With `MOSFLM_CPP_PREDICT=1` it adds `-DMOSFLM_CPP_PREDICT`, links
`libmosflm_predict_hook` + `libmosflm_geometry` via `g++` (with Homebrew
`libgfortran` path), and omits gfortran-only `FFLAGS` from the link line.

Fixed-form note: keep `CALL PREDICT_CPP_FROM_COMMON` continuations ≤72 columns
(earlier `ILG`+`XG` truncation produced a bogus `ILXG` symbol).

## COMMON mapping

| Symbol | Role | Source |
| ------ | ---- | ------ |
| `CELL(6)` | Unit cell | `/CELLCOM/` via `PREDICT_CPP_FROM_COMMON` |
| `XTOFD` | Distance, **10 µm** | `/XY/` → mm ×0.01 |
| `WAVE`, `PHIBEG` | Å, φ | `/MISC/` (in `REFLFIND`) |
| `RAST`/`RASTY` | Pixel mm | `/SCN/` (slow/fast) |
| `IYLEN`/`NREC` | Fast/slow dims | `/SCN/` |
| `DSTMAX` | → `DMINCPP=WAVE/DSTMAX` | `/REEKE/` |
| `IHG…IMG`, `XG/YG` | Spot list | `/GENDATA/` written by hook |
| `NSPOT` | Count | `/IOSP/` |

`XG`/`YG` are beam-relative **10 µm** units (mm×100). `IMG`: 0=full, 50=partial stand-in.

### Axis mapping (known-indexed)

gendata.f: `XG` ∥ lab **Y**; `YG` ∥ lab **Z**. After `mosflm_lab_to_cpp`, insulin
`XDET`/`YDET` matching (`cpp_mm_to_mosflm_xdet`, CTest `insulin_xdet`) gives
**`XG ↔ +y_mm`**, **`YG ↔ +x_mm`** (10 µm: `XG=100*y_mm`, `YG=100*x_mm`),
implemented in `predict_unpack_gendata`.

## Verification in mosflm-cpp (ProcessingLab)

```bash
ctest --test-dir mosflm-cpp/build -R 'predict_cpp_from_common|predict_cpp_hook|predict_gendata_pack'
```
