# Phase E: CCP4 drop-in for C++ IntegrationLab

Progressive drop-in of thin `bin/IntegrationLab` toward Mosflm 7.4.0
CONTROL/science. Fortran amalgamation remains `make fortran-oracle` →
`bin/ipmosflm` (not the default `make`).

| Milestone | Intent | Status |
| --------- | ------ | ------ |
| E0 | CLI / CONTROL UX (banner, RUN/GO, no-RUN warning, keyword policy) | **Done** |
| E1 | Wire AUTOINDEX into CONTROL (DPS runner) | **Done** |
| E2 | NMASKS / POSTREFINE / NADD CONTROL + `PostRefiner::refine_joint` | **Done** (joint IDXREF model still thin) |
| E3 | STRATEGY CONTROL (completeness estimate scaffold) | **Done** (oracle dump TBD) |
| E4 | MULTILATTICE store + REFIX explicit error | **Done** |
| E5 | ImageIO by extension (SMV/CBF; HDF5 message) | **Done** |
| E6 | iMosflm socket façade (`ImosflmSocketServer`) | **Done** (protocol TBD) |
| E7 | Packaging / Makefile clarity + oracle rpath helper | **Done** |

## Build

```bash
# Thin C++ (default) → IntegrationLab (.exe on Windows)
make cpp
./bin/IntegrationLab </dev/null   # banner + NO "RUN" CARD GIVEN

# Fortran oracle (does NOT use MOSFLM_FORTRAN_ORACLE=1 alone at repo root)
make dps-lib fortran-oracle       # → bin/ipmosflm
export DYLD_LIBRARY_PATH=/Applications/ccp4-9/lib
# optional: tools/fix_oracle_rpath.sh
```

## CTests

`integrationlab`, `control_session`, `control_autoindex` (if DPS), `image_io`,
`imosflm_socket`, `cbf_process_compare` (F77 vs C++ CBF PROCESS; soft-skip),
plus existing insulin PROCESS suite.
