# known_indexed image set

Real CCP4 example frame `insulin_1_001.img` (see `geometry.txt`). Scripts:

| File | Role |
| ---- | ---- |
| `cpp_process.in` | C++ `ProcessingLab` / `ControlSession` PROCESS (inline MATRIX) |
| `ipmosflm_process.in` | F77 `ipmosflm` PROCESS-only peer (no AUTOINDEX) |
| `insulin_1.mat` | NEWMAT-format orientation for F77 `MATRIX` keyword |
| `ipmosflm.in` | Full F77 oracle capture (spot search + autoindex + PROCESS) |

Performance comparison (C++ vs F77 wall time):

```bash
./tools/bench_ipmosflm_perf.sh
# MOSFLM_IPMOSFLM=/apps/ccp4/ccp4-9/bin/ipmosflm MOSFLM_PERF_REPEATS=5 ./tools/bench_ipmosflm_perf.sh
```

Oracle capture (spot list + MTZ dumps):

```bash
./tools/capture_ipmosflm_oracle.sh
```

Outputs land in `testdata/oracle/known_indexed/`.
