# Known-indexed oracle (real `ipmosflm` capture)

Captured by `tools/capture_ipmosflm_oracle.sh` from the CCP4 example frame
`insulin_1_001.img` (SMV, 2304x2304) using Mosflm 7.4.0.

| File | Contents |
| ---- | -------- |
| `ipmosflm.log` | Full run log: DPS autoindex (solution 12, I23), refined cell, orientation matrices, integration stats |
| `insulin_1_001.mtz` | Integrated frame 1; 1603 reflections. Columns: `H K L M/ISYM BATCH I SIGI IPR SIGIPR FRACTIONCALC XDET YDET ROT WIDTH LP MPART FLAG BGPKRATIOS` |
| `insulin_1.spt` | 612-space spot list from the phase-1 spot search |

Geometry / symmetry / orientation are recorded in
`../images/known_indexed/geometry.txt`.

`XDET`/`YDET` are the predicted detector **pixel** coordinates ipmosflm used for
each `H K L`; `ROT` is the reflection rotation angle. These are the raw
positional oracle for a future `PREDICT` differential test.

`predict_still.dump` / `predict_osc.dump` remain the synthetic analytical
fixtures (`tools/generate_synthetic_known_indexed.sh`, `mosflm_cpp_cli`); they
are **not** from `ipmosflm`.

## Differentials

| Test | What it checks |
| ---- | -------------- |
| `known_indexed_oracle` | Synthetic analytical dumps still match C++ fixtures |
| `insulin_hkl_overlap` | C++ PREDICT vs `insulin_1_001.mtz` after `mosflm_amat_to_cpp_ub`; ASU/Friedel-aware recall |
| `insulin_xdet` | Exact-HKL `XDET`/`YDET` residuals via `cpp_mm_to_mosflm_xdet` |

**Baseline (2026-10-17):** ASU-aware HKL recall **100%** (1463 / 1463) with
`R = [[0,0,1],[0,-1,0],[1,0,0]]`. Exact-index recall ~22% (SG/Friedel choice).
`XDET`/`YDET`: `XDET↔y_mm`, `YDET↔x_mm` relative to BEAM (94.0, 94.7) mm;
exact-HKL median residual ~**37 px** (~3 mm), inlier RMS ~38 px (no TILT/TWIST).
See `../images/known_indexed/geometry.txt`.
