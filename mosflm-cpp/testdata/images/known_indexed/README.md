# known_indexed image set (placeholder)

Place a small already-indexed dataset here, then create `ipmosflm.in` with
WAVE / DISTANCE / BEAM / CELL / MATRIX / IMAGE / GO keywords.

Then re-run:

```bash
./tools/capture_ipmosflm_oracle.sh
```

Outputs land in `testdata/oracle/known_indexed/` (log, MTZ, and — once
instrumented — a PREDICT dump matching `PREDICT_DUMP` format).
