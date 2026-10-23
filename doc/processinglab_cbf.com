#!/usr/bin/env bash
# ProcessingLab C++ CONTROL driver (--image is supported here).
./bin/ProcessingLab --image ../test_xray_data_cbf/76c1b_0002212_1_0001.cbf HKLOUT test_pl.mtz <<'EOF'
DIRECTORY ../test_xray_data_cbf
TEMPLATE 76c1b_0002212_1_####.cbf
NPROC 8
FINDSPOTS THRESHOLD 5.0 MAXSPOTS 1500 FIND 1 FILE /tmp/testdata_pl.spt MULTIFRAME RANGE 1 50
GO
AUTOINDEX DPS THRESHOLD 15
GO
PROCESS 1 TO 900
GO
EXIT
EOF
