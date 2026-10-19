# Phase C: splicing C++ PREDICT into the amalgamation

**Status: Done** (D0 closed). Insulin PROCESS with `MOSFLM_CPP_PREDICT=1` reaches
~89% exact-HKL recall; CTests `insulin_hkl_overlap` / `insulin_xdet` /
`insulin_process_mtz` gate the differentials.

## Goal

Route Fortran reflection generation through the C++ predictor (lab-frame AMAT
map + oscillation φ) without rewriting REEK.

## Hook locations (applied)

Behind `#ifdef MOSFLM_CPP_PREDICT` (default **off**):

1. **`REFLFIND`** (~194120) — interactive re-predict
2. **`GENERATE`** (MODE=0 `CALL REEK`) — PROCESS reflection list

Both call `PREDICT_CPP_FROM_COMMON`, which reads geometry from COMMON and
writes `/GENDATA/` + `/IOSP/`.

```fortran
#ifdef MOSFLM_CPP_PREDICT
        CALL PREDICT_CPP_FROM_COMMON(IFLAGPR, NREFLS,
     +       IHG, IKG, ILG, XG, YG, PHIG, PHIWG, IMG, NSPOT)
#else
        CALL PREDICT(PSI, IFLAGPR)   ! or CALL REEK(...) in GENERATE
#endif
```

## Stack

| Layer | File |
| ----- | ---- |
| Amalgamation ifdef | `mosflm/mosflm_all_ip_inc_noX.f` |
| COMMON reader | [`../fortran/predict_cpp_from_common.f90`](../fortran/predict_cpp_from_common.f90) |
| Hook | [`../fortran/predict_cpp_hook.f90`](../fortran/predict_cpp_hook.f90) |
| GENDATA pack | [`../fortran/predict_gendata_pack.f90`](../fortran/predict_gendata_pack.f90) |
| C ABI | [`../include/mosflm/legacy/MosflmCApi.h`](../include/mosflm/legacy/MosflmCApi.h) (`predict_osc_amat_cpp_`, `mosflm_ub_from_mosflm_amat`) |

## Enabling

```bash
cmake -S processing-lab -B processing-lab/build -DMOSFLM_WITH_FORTRAN=ON
cmake --build processing-lab/build --target mosflm_predict_hook
cd mosflm
rm -f mosflm_all_ip_inc_noX.o
make MOSFLM_CPP_PREDICT=1 ipmosflm
```

## COMMON mapping (read inside `PREDICT_CPP_FROM_COMMON`)

| Symbol | COMMON | Role |
| ------ | ------ | ---- |
| `CELL(6)` | `/CELLCOM/` | Unit cell |
| `XTOFD` | `/XY/` | Distance, **10 µm** → mm ×0.01 |
| `RAST`/`RASTY`, `IYLEN`/`NREC` | `/SCN/` | Pixel mm; fast/slow dims |
| `WAVE`, `PHIBEG`/`PHIEND` | `/MISC/` | Å; oscillation range |
| `AMAT`, `DSTMAX` | `/REEKE/` | Orientation [A]; `DMIN=WAVE/DSTMAX` |
| `IHG…IMG`, `XG`/`YG` | `/GENDATA/` | Written by hook |
| `NSPOT` | `/IOSP/` | Count |

Orientation: `UB_cpp = mosflm_lab_to_cpp() * AMAT`. Detector: `XG↔y_mm`,
`YG↔x_mm` (10 µm). Ewald tol fixed at 0.0025 Å⁻¹; φ step 0.05°.

## Limits vs full REEK

C++ path is fine-sliced Ewald prediction with mosaic shell width, systematic
absences, beam divergence / ice rings (REEK E0), and ``ROTMAT2`` TILT/TWIST
(``DetectorGeometry::use_rotmat2``). Insulin XDET median ~23 px with refined
fiducials; further closure needs POSTREFL orientation.

## Verification

```bash
ctest --test-dir processing-lab/build -R 'predict_cpp_from_common|predict_gendata|insulin_'
# Insulin PROCESS with CPP PREDICT:
#   IPMOSFLM=bin/ipmosflm processing-lab/tools/capture_ipmosflm_oracle.sh
```
