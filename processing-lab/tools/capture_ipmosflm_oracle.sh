#!/usr/bin/env bash
# Capture PREDICT / integrate oracles from legacy ipmosflm (ProcessingLab
# Autotools binary, or CCP4's) for a known-indexed frame.
#
# Default target is the CCP4 example frame insulin_1_001.img. Prefer
# ProcessingLab's bin/ipmosflm when built; fall back to CCP4. On Linux the
# binary may want libncurses.so.5 / libtinfo.so.5 — we synthesise .so.6
# compat symlinks under a temp dir when needed. On macOS the script sets
# DYLD_LIBRARY_PATH itself (CCP4 lib + Homebrew HDF5).
#
# Outputs land in testdata/oracle/known_indexed/ (log, MTZ, spot list) and
# testdata/images/known_indexed/ipmosflm.in is the phase-2 recipe.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IPMOSFLM_ROOT="$(cd "${ROOT}/.." && pwd)"

CCP4="${CCP4:-/Applications/ccp4-9}"
# Prefer ProcessingLab rebuild, then CCP4.
if [[ -z "${IPMOSFLM:-}" ]]; then
  if [[ -x "${IPMOSFLM_ROOT}/bin/ipmosflm" ]]; then
    IPMOSFLM="${IPMOSFLM_ROOT}/bin/ipmosflm"
  else
    IPMOSFLM="${CCP4}/bin/ipmosflm"
  fi
fi
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

export CCP4
case "$(uname -s)" in
  Darwin)
    # macOS strips DYLD_* from the environment of protected binaries such as
    # /bin/bash, so a DYLD_LIBRARY_PATH set by the caller never reaches
    # ipmosflm (it aborted in phase 1). Set it here, after bash has started:
    # CCP4's libraries, plus Homebrew HDF5 for a ProcessingLab-built binary.
    DYLD_DIRS="${CCP4}/lib"
    for hdf5 in /opt/homebrew/opt/hdf5/lib /usr/local/opt/hdf5/lib; do
      if [[ -d "${hdf5}" ]]; then
        DYLD_DIRS="${DYLD_DIRS}:${hdf5}"
        break
      fi
    done
    export DYLD_LIBRARY_PATH="${DYLD_DIRS}${DYLD_LIBRARY_PATH:+:${DYLD_LIBRARY_PATH}}"
    ;;
  *)
    # --- ncurses compat (ipmosflm wants libncurses.so.5 / libtinfo.so.5) -----
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
    export LD_LIBRARY_PATH="${COMPAT}:/usr/lib:/usr/lib/x86_64-linux-gnu:${CCP4}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
    ;;
esac

WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT
# On failure show the end of the log: the work dir is removed on exit.
fail() {
  echo "error: ipmosflm $1 failed; last lines of $2:" >&2
  tail -20 "$2" >&2 || true
  exit 1
}
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
"${IPMOSFLM}" < spots.in > spots.log 2>&1 || fail "spot search" spots.log

# --- Phase 2: autoindex from spot list + integrate one frame ----------------
cp "${SCRIPT}" process.in
echo "Phase 2: autoindex + PROCESS 1 to 1"
"${IPMOSFLM}" < process.in > process.log 2>&1 || fail "process" process.log

# --- Collect -----------------------------------------------------------------
mkdir -p "${OUT_DIR}"
cp -f process.log "${OUT_DIR}/ipmosflm.log"
[[ -f insulin_1.spt ]] && cp -f insulin_1.spt "${OUT_DIR}/insulin_1.spt"
for mtz in *.mtz; do
  [[ -f "${mtz}" ]] && cp -f "${mtz}" "${OUT_DIR}/"
done

echo "Capture finished -> ${OUT_DIR}"
echo "Inspect ipmosflm.log for the refined cell + orientation matrix."
