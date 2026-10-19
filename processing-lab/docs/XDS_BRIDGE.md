# XDS as an optional external bridge

ipmosflm and mosflm-cpp **do not** require [XDS](https://wiki.uni-konstanz.de/xds/index.php/Installation) to build or run. XDS is useful when you need to **convert** an unsupported detector format into miniCBF (or process HDF5 with vendor plugins), then feed the result into MOSFLM’s existing CBFlib path ([CBF_STATUS.md](CBF_STATUS.md)).

## When to use it

| Situation | Suggested workflow |
| --------- | ------------------ |
| Native CBF / miniCBF | Read directly (`DETECTOR CBF` / Pilatus / `cbfwrap_`) |
| Format XDS understands but MOSFLM does not | Convert with XDS `2cbf` (or related tools) → miniCBF → MOSFLM |
| Eiger HDF5 | Prefer MOSFLM `fromhdf5`, **or** convert frames externally (`H5ToXds` / `hdf2mini-cbf`) if needed |
| `*.cbf.gz` | XDS `xds-zcbf` plugin can speed XDS; for MOSFLM, decompress or convert first |

## Install notes (operators only)

Follow the [XDSwiki Installation](https://wiki.uni-konstanz.de/xds/index.php/Installation) page:

- Academic download via Dectris form; commercial users need a subscription.
- macOS Apple Silicon: `XDS-Apple_M1.zip` → `/usr/local/bin` or `$HOME/bin`.
- Quarantine: `sudo xattr -dr com.apple.quarantine XDS-Apple_M1/*` if Gatekeeper blocks binaries.
- HDF5 plugins (Durin / Neggia) and helpers (`generate_XDS.INP`, XDSGUI) are optional extras from the wiki `get_folder_silicon.sh` scripts.
- CCP4 helps XDSSTAT / XDSCONV / POINTLESS from XDSGUI; it is **not** required to compile ipmosflm’s CBF reader.

## Explicit non-goals

- Do **not** vendor XDS binaries into this repository.
- Do **not** add Neggia/Durin as link-time dependencies of `bin/ProcessingLab`,
  `bin/ipmosflm` (oracle), or `processing-lab`.
- Do **not** treat `2cbf` as a runtime requirement for every image open.

## Importing Xia2 `.xinfo` / `XDS.INP`

ProcessingLab can seed CONTROL from Xia2/XDS text files without running XDS:

```bash
./ProcessingLab --xinfo automatic.xinfo --xds-inp path/to/XDS.INP HKLOUT out.mtz script.in
```

| Flag | Source | Maps into CONTROL |
|------|--------|-------------------|
| `--xds-inp` | `XDS.INP` | `WAVELENGTH`, `DISTANCE`, `PIXEL`, `BEAM` (from `ORGX`/`ORGY` x `QX`/`QY` when present), `CELL`, `SYMMETRY`, `DIRECTORY`+`TEMPLATE` (`?` becomes `#`), `PROCESS` (+ `ANGLE` from `OSCILLATION_RANGE`) |
| `--xinfo` | Xia2 `.xinfo` | `WAVELENGTH`, `DIRECTORY`, `IMAGE`/`TEMPLATE`, `PROCESS` from first sweep |

Apply order is **XDS.INP then `.xinfo`**, so local Xia2 paths and `START_END` override stale absolute `NAME_TEMPLATE_OF_DATA_FRAMES` paths.

Orientation is not in `XDS.INP` / `.xinfo`; use `--xparm XPARM.XDS` (or CONTROL `XPARM` / `GXPARM`) to seed `CELL`, `SYMMETRY`, `DISTANCE`, `PIXEL`, `BEAM`, and `MATRIX` when present. `AUTOINDEX` GO writes `XPARM.XDS`; `PROCESS` GO writes `GXPARM.XDS` (disable with `XPARMOUT OFF`). Prefer an index-stage `XDS.INP` (or the image header) for beam centre when the integrate/CORRECT `XDS.INP` omits `ORGX`/`ORGY`.

CORRECT-lite: `POLARIZATION ρ` (default 0.98) in LP; `AIR` path-length absorption; systematic-absence I/σ sanity; `CELLFIX` before PROCESS.

INTEGRATE-lite: `PARTIALSCALE` (I÷FRACTIONCALC), `AUTOSIGMA` (σM/σD from partials), PEAKHALF clamped to half SEPARATION, trusted-mask profile bg.

### IDXREF / AUTOINDEX (COLSPOT coverage)

XDS IDXREF uses multi-image COLSPOT (often thousands of spots). ProcessingLab
`FINDSPOTS` defaults to `MAXSPOTS 1500`. For large cells, **single-frame
`MAXSPOTS ≤ 500` is unsafe** — DPS can lock onto a short-vector / sublattice
(e.g. 76c1b wrongly indexing as ~43×118×60 mC instead of ~62×187×72 mP).

Recommended CONTROL for indexing:

```text
FINDSPOTS THRESHOLD 5 MAXSPOTS 1500 FIND 1 FILE spots.spt MULTIFRAME RANGE 1 50
GO
AUTOINDEX DPS THRESHOLD 15 …
GO
```

`MULTIFRAME` without `RANGE` defaults to a 50-frame wedge from `FIND`.
For CBF/SMV templates, each frame path is resolved separately (HDF5 masters
still use the frame index). Default `MERGE_XY 0` keeps per-frame 3D spots;
set `MERGE_XY` / `MERGE_PHI` only when clustering persistent peaks.
DiffractionLab `build_findspots_commands(..., for_indexing=True)` emits the
same MULTIFRAME pattern with `MAXSPOTS 1500`.

After DPS + Bravais SUGGLATT, ProcessingLab applies a BESTSOLAR-lite gate
(`%indexed` / RMSD HKL) and may fall back among LASTGOOD candidates (preferring
`mP` over a poorly indexing `mC`). Look for `INDEXFIT:` / `PROCESSINGLAB_INDEXED`
in the log. If AUTOINDEX still fails, seed from Xia2/XDS with `--xparm`
pointing at an **index-stage** `XPARM.XDS` (not a prior bad ProcessingLab
`GXPARM.XDS`).

Optional multi-image COLSPOT: `FINDSPOTS MULTIFRAME RANGE first last …`. Multi-frame PROCESS runs periodic post-refinement (`POSTREFINE PERIOD n`, default every 25 frames).
