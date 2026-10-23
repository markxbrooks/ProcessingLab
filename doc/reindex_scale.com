
#!/usr/bin/env bash
set -euo pipefail

FILESTEM="process"
PROG1="pointless"
PROG2="aimless"

# Determine input and output filenames
POINTLESS_MTZ="${FILESTEM}_${PROG1}.mtz"
AIMLESS_MTZ="${FILESTEM}_${PROG2}.mtz"

# Determine the space group and reindex if necessary
"${PROG1}" \
    HKLIN "${FILESTEM}.mtz" \
    HKLOUT "${POINTLESS_MTZ}"

# Scale and merge reflections
"${PROG2}" \
    HKLIN "${POINTLESS_MTZ}" \
    HKLOUT "${AIMLESS_MTZ}" << EOF
GO
EOF
