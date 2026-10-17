#!/usr/bin/env bash
# Capture PREDICT / integrate oracles from ipmosflm for a known-indexed frame.
#
# Default target is the CCP4 example frame insulin_1_001.img. ipmosflm is a
# Fortran binary that links libncurses.so.5 / libtinfo.so.5; on systems that
# only ship .so.6 we synthesise compatibility symlinks under a temp dir.
#
# Outputs land in testdata/oracle/known_indexed/ (log, MTZ, spot list) and
# testdata/images/known_indexed/ipmosflm.in is the phase-2 recipe.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IPMOSFLM_ROOT="$(cd "${ROOT}/.." && pwd)"

CCP4="${CCP4:-/apps/ccp4/ccp4-9}"
IPMOSFLM="${IPMOSFLM:-${CCP4}/bin/ipmosflm}"
IMAGE="${MOSFLM_IMAGE:-${CCP4}/examples/data/insulin_1_001.img}"
TEMPLATE="${MOSFLM_TEMPLATE:-insulin_1_###.img}"

OUT_DIR="${ROOT}/testdata/oracle/known_indexed"
IMG_DIR="${ROOT}/testdata/images/known_indexed"
SCRIPT="${IMG_DIR}/ipmosflm.in"

if [[ ! -x "${IPMOSFLM}" ]]; then
  echo "error: ipmosflm not found at ${IPMOSFLM}" >&2
  exit 1
fi
if [[ ! -f "${IMAGE}" ]]; then
  echo "error: input image not found at ${IMAGE}" >&2
  echo "set MOSFLM_IMAGE to a known-indexed frame" >&2
  exit 1
fi

# --- ncurses compat (ipmosflm wants libncurses.so.5 / libtinfo.so.5) ---------
COMPAT="${TMPDIR:-/tmp}/opencode_ccp4compat"
mkdir -p "${COMPAT}"
for lib in ncurses tinfo; do
  if [[ ! -e "${COMPAT}/lib${lib}.so.5" ]]; then
    for cand in /usr/lib/x86_64-linux-gnu "/usr/lib64" "/lib/x86_64-linux-gnu"; do
      if [[ -e "${cand}/lib${lib}.so.6" ]]; then
        ln -sf "${cand}/lib${lib}.so.6" "${COMPAT}/lib${lib}.so.5"
        break
      fi
    done
  fi
done

export CCP4
export LD_LIBRARY_PATH="${COMPAT}:/usr/lib:/usr/lib/x86_64-linux-gnu${CCP4}/lib:${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT
cd "${WORK}"
ln -sf "${IMAGE}" "$(basename "${IMAGE}")"

# --- Phase 1: spot finding (writes insulin_1.spt) ---------------------------
cat > spots.in <<EOF
IMAGE $(basename "${IMAGE}")
DISTANCE 160.0
WAVELENGTH 0.9790
PIXEL 0.0816 0.0816
BEAM 94.0 94.7
SEPARATION 0.57 0.57
AUTOINDEX
GO
EOF
echo "Phase 1: spot search -> insulin_1.spt"
"${IPMOSFLM}" < spots.in > spots.log 2>&1 || {
  echo "ipmosflm spot search failed; see ${WORK}/spots.log" >&2; exit 1; }

# --- Phase 2: autoindex from spot list + integrate one frame ----------------
cp "${SCRIPT}" process.in
echo "Phase 2: autoindex + PROCESS 1 to 1"
"${IPMOSFLM}" < process.in > process.log 2>&1 || {
  echo "ipmosflm process failed; see ${WORK}/process.log" >&2; exit 1; }

# --- Collect -----------------------------------------------------------------
mkdir -p "${OUT_DIR}"
cp -f process.log "${OUT_DIR}/ipmosflm.log"
[[ -f insulin_1.spt ]] && cp -f insulin_1.spt "${OUT_DIR}/insulin_1.spt"
for mtz in *.mtz; do
  [[ -f "${mtz}" ]] && cp -f "${mtz}" "${OUT_DIR}/"
done

echo "Capture finished -> ${OUT_DIR}"
echo "Inspect ipmosflm.log for the refined cell + orientation matrix."
