# Phase E: CCP4 drop-in for C++ ProcessingLab

Progressive drop-in of thin `bin/ProcessingLab` toward Mosflm 7.4.0
CONTROL/science. Fortran amalgamation remains `make fortran-oracle` →
`bin/ipmosflm` (not the default `make`).

| Milestone | Intent | Status |
| --------- | ------ | ------ |
| E0 | CLI / CONTROL UX (banner, RUN/GO, no-RUN warning, keyword policy) | **Done** |
| E1 | Wire AUTOINDEX into CONTROL (DPS runner) | **Done** |
| E2 | NMASKS / POSTREFINE / NADD CONTROL + `PostRefiner::refine_joint` | **Done** (joint IDXREF model still thin) |
| E3 | STRATEGY CONTROL (completeness estimate scaffold) | **Done** (oracle dump TBD) |
| E4 | MULTILATTICE store + REFIX explicit error | **Done** |
| E5 | ImageIO by extension (SMV/CBF; HDF5 decode for Eiger/DLS masters) | **Done** |
| E6 | iMosflm socket façade (`ImosflmSocketServer`) | **Done** (protocol TBD) |
| E7 | Packaging / Makefile clarity + oracle rpath helper | **Done** |
| E8 | PROCESS parity vs F77 (REEK-ish predict, 18-col MTZ, ASU/`M/ISYM`, `FLAG`) | **Done** (soft count gate; tighten toward 1603 rows) |

## Build

```bash
# Thin C++ (default) → ProcessingLab (.exe on Windows)
make cpp
./bin/ProcessingLab </dev/null   # banner + NO "RUN" CARD GIVEN

# Fortran oracle (does NOT use MOSFLM_FORTRAN_ORACLE=1 alone at repo root)
make dps-lib fortran-oracle       # → bin/ipmosflm
export DYLD_LIBRARY_PATH=/Applications/ccp4-9/lib
# optional: tools/fix_oracle_rpath.sh
```

## CTests

`integrationlab`, `control_session`, `control_autoindex` (if DPS), `image_io`,
`imosflm_socket`, `cbf_process_compare` (F77 vs C++ CBF PROCESS; soft-skip),
`smv_process_compare` (F77 vs ProcessingLab insulin PROCESS — cell/SG, ASU +
exact HKL recall, soft count/I-ratio; release gate for E8),
`smv_dials_compare` (DIALS vs ProcessingLab on the same SMV — cell/SG, XDET
proximity, soft count/I-ratio; soft-skip without DIALS), `image_io` (SMV +
optional `MOSFLM_HDF5_FIXTURE` for Eiger masters), plus existing insulin
PROCESS suite (`insulin_cpp_process`, `insulin_control`, `insulin_isi`,
`insulin_xdet`, `insulin_profile`).

### HDF5 / NeXus masters

With HDF5 at configure time, ProcessingLab can load Dectris Eiger
`*_master.h5` (and DLS-style `.nxs`) via `Hdf5ImageReader`. CONTROL
`PROCESS` first-image selects the 1-based frame. External `data_*.h5`
chunks must sit next to the master. Example (script still required for
`GO`):

```bash
./bin/ProcessingLab HKLOUT /tmp/out.mtz \
  --image /path/to/thaumatin_9_1_master.h5 your_process.in
```

### AUTOINDEX Bravais (SUGGLATT / LASTGOOD)

After DPS, ProcessingLab selects a lattice like Mosflm `SUGGLATT` /
`DPSCHOOSE` (penalty-only):

1. Build Mosflm MERIT from the raw DPS goof for jump math only
   (`round(crit * 999 / max_crit)`).
2. Find `LASTGOOD` — the sorted rank before the largest relative MERIT jump
   among candidates with MERIT ≤ 100 (same idea as Mosflm `TESTGOOD`).
3. Among ranks `1..LASTGOOD` (when a jump was found), prefer the highest
   crystal class (`c` > `h`/`r` > `t` > `o` > `m` > `a`).
4. If no jump is detected, fall back to the soft cut (criterion ≤ 50) with
   the same crystal-class ranking.

Then apply that Bravais OM to CELL/MATRIX. Default `SYMMETRY` follows Mosflm
`SPGMIN` (e.g. `tP` → P4 / #75) when unset. Overrides: `AUTOINDEX … SYSTEM tP`
or `prefer_character`. (Full Mosflm prerefine/RMSD `BESTSOLAR` is not ported.)

AUTOINDEX also prints the Mosflm-style lattice-character table
(`No PENALTY LATT a b c alpha beta gamma Possible spacegroups`), sorted by
penalty. PENALTY is the raw DPS Niggli goof (not MERIT 0–999). The table
notes `LASTGOOD` when a jump is found and ends with a `Suggested Solution`
footer.

### FINDSPOTS → `.spt` → AUTOINDEX

ProcessingLab includes a DiffractionLab-style peak finder (`SpotFinder`:
MAD background + 3×3 local maxima) and MOSFLM `.spt` write/read. CONTROL:

```text
FINDSPOTS THRESHOLD 5.0 MAXSPOTS 500 FIND 1 PHI 0.0 0.15 FILE thaumatin.spt
GO
AUTOINDEX DPS THRESHOLD 20 FILE thaumatin.spt
GO
```

`--image` supplies the frame (SMV/CBF/HDF5). If `AUTOINDEX` omits `FILE`,
the path from the preceding `FINDSPOTS` is reused. Geometry cards may be
omitted when the image header fills distance / wavelength / beam / pixel.

### Multi-frame PROCESS

`PROCESS 1 TO N` / `PROCESS 1 N` loops images on one `GO`. HDF5/NeXus masters
reuse the path with a 1-based frame index; SMV/CBF series need `TEMPLATE`.
Each frame’s reflections get `BATCH =` image number; φ comes from the image
header unless `PROCESS … START` / `ANGLE` is set. Example:

```bash
./bin/ProcessingLab HKLOUT /tmp/out.mtz \
  --image /path/to/master.h5 <<'EOF'
FINDSPOTS THRESHOLD 5 FIND 1 FILE /tmp/spots.spt
GO
AUTOINDEX DPS THRESHOLD 20
GO
PROCESS 1 TO 10
GO
EXIT
EOF
```

## DiffractionLab

Sibling app [`DiffractionLab`](../DiffractionLab) can drive `ProcessingLab` via
CONTROL scripts (`diffraction_lab/processinglab_client.py`). After AUTOINDEX,
`ProcessingLab` prints `PROCESSINGLAB_CELL` / `PROCESSINGLAB_MATRIX` (and
writes `integrationlab_index.txt`). Set `DIFFRACTIONLAB_BACKEND=processinglab`
or `DIFFRACTIONLAB_PREFER_PROCESSINGLAB=1`.
