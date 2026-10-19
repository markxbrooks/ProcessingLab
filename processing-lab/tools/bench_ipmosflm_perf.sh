#!/usr/bin/env bash
# Wall-time comparison: ProcessingLab (C++ / ProcessingLab) vs F77 ipmosflm.
#
# Defaults:
#   C++  → $MOSFLM_PROCESSINGLAB or build/ProcessingLab / ../bin/ProcessingLab
#   F77  → $MOSFLM_IPMOSFLM or /apps/ccp4/ccp4-9/bin/ipmosflm ($CCP4/bin/ipmosflm)
#   IMG  → $MOSFLM_IMAGE or $CCP4/examples/data/insulin_1_001.img
#
# Usage:
#   ./tools/bench_ipmosflm_perf.sh
#   MOSFLM_PERF_REPEATS=5 ./tools/bench_ipmosflm_perf.sh
#   MOSFLM_IPMOSFLM=/apps/ccp4/ccp4-9/bin/ipmosflm \
#     MOSFLM_PROCESSINGLAB=./build/ProcessingLab ./tools/bench_ipmosflm_perf.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${ROOT}"

CCP4="${CCP4:-/apps/ccp4/ccp4-9}"
export CCP4

if [[ -z "${MOSFLM_PROCESSINGLAB:-}" ]]; then
  for cand in \
    "${ROOT}/build/ProcessingLab" \
    "${ROOT}/../bin/ProcessingLab" \
    "${ROOT}/bin/ProcessingLab"
  do
    if [[ -x "${cand}" ]]; then
      MOSFLM_PROCESSINGLAB="${cand}"
      break
    fi
  done
fi

if [[ -z "${MOSFLM_IPMOSFLM:-}" ]]; then
  if [[ -x /apps/ccp4/ccp4-9/bin/ipmosflm ]]; then
    MOSFLM_IPMOSFLM=/apps/ccp4/ccp4-9/bin/ipmosflm
  elif [[ -x "${CCP4}/bin/ipmosflm" ]]; then
    MOSFLM_IPMOSFLM="${CCP4}/bin/ipmosflm"
  elif [[ -x "${ROOT}/../bin/ipmosflm.oracle" ]]; then
    MOSFLM_IPMOSFLM="${ROOT}/../bin/ipmosflm.oracle"
  elif [[ -x "${ROOT}/../bin/ipmosflm" ]]; then
    MOSFLM_IPMOSFLM="${ROOT}/../bin/ipmosflm"
  fi
fi

export MOSFLM_PROCESSINGLAB MOSFLM_IPMOSFLM
export MOSFLM_IMAGE="${MOSFLM_IMAGE:-${CCP4}/examples/data/insulin_1_001.img}"
export MOSFLM_PERF_REPEATS="${MOSFLM_PERF_REPEATS:-3}"

BIN=""
for cand in \
  "${ROOT}/build/test_ipmosflm_perf" \
  "${ROOT}/../cmake-build-debug/mosflm-cpp/test_ipmosflm_perf"
do
  if [[ -x "${cand}" ]]; then
    BIN="${cand}"
    break
  fi
done

if [[ -z "${BIN}" ]]; then
  echo "building test_ipmosflm_perf (and ProcessingLab if needed)..."
  cmake -S "${ROOT}" -B "${ROOT}/build" -DCMAKE_BUILD_TYPE=Release
  cmake --build "${ROOT}/build" --target test_ipmosflm_perf ProcessingLab -j"$(nproc 2>/dev/null || echo 4)"
  BIN="${ROOT}/build/test_ipmosflm_perf"
  export MOSFLM_PROCESSINGLAB="${MOSFLM_PROCESSINGLAB:-${ROOT}/build/ProcessingLab}"
fi

echo "C++ binary : ${MOSFLM_PROCESSINGLAB:-unset}"
echo "F77 binary : ${MOSFLM_IPMOSFLM:-unset}"
echo "Image      : ${MOSFLM_IMAGE}"
echo "Repeats    : ${MOSFLM_PERF_REPEATS}"
echo

exec "${BIN}"
