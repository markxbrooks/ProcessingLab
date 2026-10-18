#!/usr/bin/env bash
# Point a locally built Fortran oracle ipmosflm at CCP4-9 dylibs so it runs
# without DYLD_LIBRARY_PATH (Phase E / E7). macOS only.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="${1:-$ROOT/bin/ipmosflm}"
CCP4_LIB="${CCP4:-/Applications/ccp4-9}/lib"

if [[ ! -x "$BIN" ]]; then
  echo "fix_oracle_rpath: missing executable $BIN" >&2
  exit 1
fi
if [[ ! -d "$CCP4_LIB" ]]; then
  echo "fix_oracle_rpath: missing CCP4 lib dir $CCP4_LIB" >&2
  exit 1
fi

# Only rewrite if this looks like the amalgamation (~MBs), not thin C++.
size=$(stat -f%z "$BIN" 2>/dev/null || stat -c%s "$BIN")
if (( size < 2000000 )); then
  echo "fix_oracle_rpath: $BIN looks like thin C++ ($size bytes); refusing" >&2
  exit 2
fi

for lib in libccp4f.8.0.dylib libccp4c.8.0.dylib libmmdb2.0.dylib; do
  if [[ -f "$CCP4_LIB/$lib" ]]; then
    # Replace any install name whose basename matches.
    old=$(otool -L "$BIN" | awk '/'"$lib"'/ {print $1; exit}')
    if [[ -n "${old:-}" ]]; then
      echo "  $old → $CCP4_LIB/$lib"
      install_name_tool -change "$old" "$CCP4_LIB/$lib" "$BIN"
    fi
  fi
done

echo "Done. Try: $BIN </dev/null | head"
