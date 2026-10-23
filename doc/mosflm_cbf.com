#!/usr/bin/env bash
# F77 ipmosflm (CCP4). Does NOT support ProcessingLab's --image switch.
#
# Do not use IMAGE here: Mosflm forbids IMAGE + PROCESS in the same run.
# TEMPLATE + FINDSPOTS FIND n is enough; PROCESS then reuses the same sweep
# and the cell/matrix from AUTOINDEX.
#
# Extra EXIT lines: after PROCESS Mosflm prints a long summary then reads
# stdin again; a single EXIT in a heredoc often hits "READ after EOF".
set -euo pipefail
SCRIPT="$(mktemp -t mosflm_cbf.XXXXXX)"
trap 'rm -f "$SCRIPT"' EXIT
cat >"$SCRIPT" <<'EOF'
DIRECTORY ../test_xray_data_cbf
TEMPLATE 76c1b_0002212_1_####.cbf
FINDSPOTS THRESHOLD 5.0 MAXSPOTS 300 FIND 1 PHI 0.0 0.2 FILE /tmp/testdata_mosflm.spt
GO
AUTOINDEX DPS THRESHOLD 20 FILE /tmp/testdata_mosflm.spt
GO
PROCESS 1 TO 900
GO
EXIT
EXIT
EOF
ipmosflm HKLOUT test_mosflm.mtz <"$SCRIPT"
