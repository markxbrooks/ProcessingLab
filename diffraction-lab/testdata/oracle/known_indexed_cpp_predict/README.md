# Insulin PROCESS with `MOSFLM_CPP_PREDICT=1`

Captured 2026-10-17 from ProcessingLab `bin/ipmosflm` built with
`make MOSFLM_CPP_PREDICT=1` (GENERATE → `PREDICT_CPP_FROM_COMMON` instead of
REEK). Same `ipmosflm.in` recipe as `../known_indexed/` (solution 12, I23).

| File | Notes |
| ---- | ----- |
| `ipmosflm_cpp_predict.log` | Full log |
| `insulin_1_001.mtz` | 4819 reflections (oracle Fortran REEK: 1603) |
| `insulin_1.spt` | Spot list copied from Fortran oracle capture |

## Vs Fortran oracle (`../known_indexed/`)

| Metric | Value |
| ------ | ----- |
| Oracle unique HKL | 1463 |
| CPP unique HKL | 4231 |
| Exact overlap / recall | 1308 / **89.4%** |

CPP path predicts more spots (no REEK mosaicity/divergence gating). Overlaps are
exact-index (same ASU choice as integrated MTZ for this run).
