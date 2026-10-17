#!/usr/bin/env bash
# Scaffold for capturing PREDICT / integrate oracles from bin/ipmosflm.
#
# Requires a known-indexed image set (CELL + MATRIX + images). Until that
# dataset is under testdata/images/known_indexed/, this script exits with a
# clear message. See testdata/README.md.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IPMOSFLM_ROOT="$(cd "${ROOT}/.." && pwd)"
IPMOSFLM="${IPMOSFLM:-${IPMOSFLM_ROOT}/bin/ipmosflm}"
OUT_DIR="${ROOT}/testdata/oracle/known_indexed"
IMG_DIR="${ROOT}/testdata/images/known_indexed"
SCRIPT="${IMG_DIR}/ipmosflm.in"

export CCP4="${CCP4:-/Applications/ccp4-9}"
export DYLD_LIBRARY_PATH="${CCP4}/lib:/opt/homebrew/opt/hdf5/lib${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"

if [[ ! -x "${IPMOSFLM}" ]]; then
  echo "error: ipmosflm not found at ${IPMOSFLM}" >&2
  exit 1
fi

if [[ ! -f "${SCRIPT}" ]]; then
  mkdir -p "${IMG_DIR}" "${OUT_DIR}"
  cat > "${IMG_DIR}/README.md" <<'EOF'
# known_indexed image set (placeholder)

Place a small already-indexed dataset here, then create `ipmosflm.in` with
WAVE / DISTANCE / BEAM / CELL / MATRIX / IMAGE / GO keywords.

Then re-run:

```bash
./tools/capture_ipmosflm_oracle.sh
```

Outputs land in `testdata/oracle/known_indexed/` (log, MTZ, and — once
instrumented — a PREDICT dump matching `PREDICT_DUMP` format).
EOF
  echo "No ${SCRIPT} yet — wrote ${IMG_DIR}/README.md scaffold and exiting."
  echo "Analytical oracles remain: testdata/oracle/predict_still.dump, predict_osc.dump"
  exit 0
fi

mkdir -p "${OUT_DIR}"
LOG="${OUT_DIR}/ipmosflm.log"
echo "Running ${IPMOSFLM} with ${SCRIPT} → ${LOG}"
"${IPMOSFLM}" < "${SCRIPT}" > "${LOG}" 2>&1 || {
  echo "ipmosflm failed; see ${LOG}" >&2
  exit 1
}
echo "Capture finished. Inspect ${OUT_DIR} and convert any prediction lists"
echo "into PREDICT_DUMP format for CTest (see docs/FIRST_ROUTINE.md)."
