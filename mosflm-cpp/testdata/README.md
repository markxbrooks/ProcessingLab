# Regression corpus for mosflm-cpp

## Layout

| Path | Purpose |
| ---- | ------- |
| `analytical/` | Geometry fixtures that need no images (INVCEL / METRIC / UB) |
| `oracle/` | Outputs from `bin/ipmosflm` on the known-indexed image set |

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

```bash
export CCP4=/Applications/ccp4-9
export DYLD_LIBRARY_PATH="$CCP4/lib:/opt/homebrew/opt/hdf5/lib"
IPMOSFLM=../../bin/ipmosflm

# Example skeleton (edit paths / keywords for your images):
$IPMOSFLM <<EOF > testdata/oracle/known_indexed/ipmosflm.log
! WAVE … DISTANCE … BEAM … CELL … MATRIX …
! GO
EXIT
EOF
```

Compare C++ outputs with **numeric tolerances**, not bit-identical floats.

## Analytical fixture

See `analytical/thaumatin_cell.txt` — used by `test_metric_invcel`.
