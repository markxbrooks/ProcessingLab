# Oracle outputs

## `predict_still.dump`

Analytical **stills** PREDICT regression dump (Ewald-sphere reference).

| Item | Value |
| ---- | ----- |
| Generator | [`tools/predict_still_oracle.py`](../../tools/predict_still_oracle.py) |
| Inputs | [`predict_still.inputs.txt`](predict_still.inputs.txt) |
| Format | `PREDICT_DUMP n=N` then `h k l x_mm y_mm phi_beg phi_end full` |

Regenerate:

```bash
python3 tools/predict_still_oracle.py > testdata/oracle/predict_still.dump
```

**Provenance note:** A full `bin/ipmosflm` PREDICT capture needs a known-indexed
image set under `testdata/images/` (see parent `testdata/README.md`). Until that
dataset is available, this independent Python Ewald stills dump is the frozen
oracle for `ReflectionPredictor` differentials. Geometry constants match the
C++ `UnitCell` / `DetectorGeometry` conventions.

## `predict_osc.dump`

Analytical **oscillation** PREDICT dump (`φ ∈ [0, 0.5]°`, `dmin = 5 Å`).

| Item | Value |
| ---- | ----- |
| Generator | [`tools/predict_osc_oracle.py`](../../tools/predict_osc_oracle.py) |
| Classification | `full=0` if on-sphere at either oscillation endpoint (partial) |

```bash
python3 tools/predict_osc_oracle.py > testdata/oracle/predict_osc.dump
```

## `integrate_synthetic.dump`

Synthetic summation-integrate intensities for two planted flat boxes (bg=0):

| hkl | peak/px | npix | I |
| --- | ------- | ---- | - |
| 1 0 0 | 50 | 25 | 1250 |
| 0 1 0 | 80 | 25 | 2000 |

Exercised by CTest `summation_integrate` (also MTZ write/read compare).

## `ipmosflm_smoke.log`

Smoke run of `bin/ipmosflm` with immediate `EXIT` (build sanity only).
