# Regression corpus for processing-lab

## Layout

| Path | Purpose |
| ---- | ------- |
| `analytical/` | Geometry fixtures that need no images (INVCEL / METRIC / UB) |
| `oracle/` | Frozen dumps for differentials (`predict_still.dump` analytical stills; future `ipmosflm` captures) |

## Known-indexed image set (to capture)

Place a small real dataset here (or symlink) once available:

1. **Simple known cell** — e.g. thaumatin-like tetragonal (`58 58 150 90 90 90`) with a fixed orientation matrix.
2. **Several consecutive images** covering partial reflections (0.5–1.0° oscillation).
3. **Weak / high-background** case for integrator stress.
4. Optional: alternate detector header (CBF or SMV) for I/O coverage.

Recommended location once obtained:

```
testdata/images/known_indexed/
  README.md          # wavelength, distance, beam centre, CELL, MATRIX
  image_001.* …      # few frames
testdata/oracle/known_indexed/
  ipmosflm.log
  integrate.mtz
  refined_params.txt
```

## Capture oracle results

Synthetic corpus (analytical C++ dumps, no real images yet):

```bash
./tools/generate_synthetic_known_indexed.sh
```

Prefer the scaffold script for a future real `ipmosflm` capture (creates
`images/known_indexed/` README if missing):

```bash
./tools/capture_ipmosflm_oracle.sh
```

Or run manually once `testdata/images/known_indexed/ipmosflm.in` exists:

```bash
export CCP4=/Applications/ccp4-9
export DYLD_LIBRARY_PATH="$CCP4/lib:/opt/homebrew/opt/hdf5/lib"
../../bin/ipmosflm < testdata/images/known_indexed/ipmosflm.in \
  > testdata/oracle/known_indexed/ipmosflm.log
```

Compare C++ outputs with **numeric tolerances**, not bit-identical floats.
Analytical PREDICT dumps (`oracle/predict_still.dump`, `predict_osc.dump`) are
used until a known-indexed capture is available.

## Analytical fixture

See `analytical/thaumatin_cell.txt` — used by `test_metric_invcel`.
