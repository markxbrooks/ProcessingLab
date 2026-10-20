# ProcessingLab

Progressive C++ port of the MOSFLM diffraction-processing core.

![ProcessingLab Logo](resources/ProcessingLabLogo.png)


**Phase D Done:** default build/install is the C++ core. The release binary is
**`ProcessingLab`** (`ProcessingLab.exe` on Windows) — a thin CONTROL →
PROCESS driver that does **not** link the Fortran amalgamation. Autotools
`mosflm/` remains an optional scientific oracle named **`ipmosflm`** (F77).

| Path | Role |
| ---- | ---- |
| [`processing-lab`](processing-lab/) | C++ geometry, PREDICT, integration, CONTROL, CTest |
| [`mosflm/`](mosflm/) | Amalgamated Fortran (**opt-in oracle** → `bin/ipmosflm`) |
| [`porting.md`](porting.md) | Migration roadmap and current status |
| `bin/ProcessingLab` | C++ CONTROL binary (default; `.exe` on Windows) |
| `bin/mosflm_cpp_cli` | Full C++ CLI |
| `bin/ipmosflm` | F77 oracle only (`make fortran-oracle`) |

## Quick start (default)

```bash
make cpp          # builds cbf/lib/libcbf.a then processing-lab
make test         # ctest in processing-lab/build
./bin/ProcessingLab version
./bin/mosflm_cpp_cli process-insulin
./bin/ProcessingLab HKLOUT /tmp/out.mtz --image /path/to/insulin_1_001.img \
  processing-lab/testdata/images/known_indexed/cpp_process.in
# C++ vs F77 PROCESS wall-time (CCP4 ipmosflm by default):
./processing-lab/tools/bench_ipmosflm_perf.sh
```

If configuring `processing-lab` directly, build CBFlib first: `make -C cbf`.

Top-level CMake (`-S .`) adds `processing-lab` when `PROCESSINGLAB_BUILD_CPP=ON`
(default).

## Optional Fortran oracle

Requires CCP4. See [`porting.md`](porting.md) for macOS notes.

```bash
make fortran-oracle   # → bin/ipmosflm (F77; does not overwrite ProcessingLab)
# or:
#   ./configure --with-ccp4-prefix=/Applications/ccp4-9
#   make -C mosflm MOSFLM_FORTRAN_ORACLE=1 ipmosflm
# Experimental C++ PREDICT splice into the oracle:
#   cmake -S processing-lab -B processing-lab/build -DMOSFLM_FORTRAN_ORACLE=ON
#   cmake --build processing-lab/build --target mosflm_predict_hook
#   make -C mosflm MOSFLM_FORTRAN_ORACLE=1 MOSFLM_CPP_PREDICT=1 ipmosflm
```

# Usage

```commandline
./bin/ProcessingLab --image ~/Downloads/thaumatin_9_1_master.h5 <<'EOF'
NPROC 1
DIRECTORY /tmp
FINDSPOTS THRESHOLD 5.0 MAXSPOTS 300 FIND 1 PHI 0.0 0.15 FILE /tmp/thaumatin_cpp.spt
GO
AUTOINDEX DPS THRESHOLD 20
GO
PROCESS 1 TO 999
GO
EXIT
EOF
```

```commandline
./bin/ProcessingLab --image ../test_xray_data_cbf/76c1b_0002212_1_0001.cbf hklout test.mtz <<'EOF'
DIRECTORY test_xray_data_cbf
TEMPLATE 76c1b_0002212_1_####.cbf
DMAX 320
FINDSPOTS THRESHOLD 5.0 MAXSPOTS 1500 FIND 1 FILE /tmp/testdata.spt MULTIFRAME RANGE 1 50
GO
AUTOINDEX DPS THRESHOLD 15 SEPMIN 6 CLUSTER_RADIUS 2
GO
PROCESS 1 TO 900
GO
EXIT
EOF
```

## License

BSD-style — see [`COPYING`](COPYING) (A.G.W. Leslie / MRC LMB).
