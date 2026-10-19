#!/bin/tcsh -f

set tttt = $argv[1]
set output = ${tttt:r}.sca

xds2sca << EOF
XDSFILE $tttt
SCAFILE $output
EOF
