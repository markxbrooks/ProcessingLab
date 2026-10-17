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
| `insulin_hkl_overlap` | C++ PREDICT vs `insulin_1_001.mtz` Miller-index set, probing A-matrix axis/sign permutations |

**Baseline (2026-10-17):** best recall ≈ **6.8%** (99 / 1463 unique MTZ HKLs)
with axis permutation **YZX** and sign pattern **2**, dmin=2.5 Å, φ∈[0,1]°.
Predictor still emits many extra spots (SG absences / Ewald window). Positional
`XDET`/`YDET` comparison remains blocked on lab-frame reconciliation — see
`../images/known_indexed/geometry.txt`.
