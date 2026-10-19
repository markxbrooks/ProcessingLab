# CBF support status (ipmosflm)

Classification relative to [critique.md](../../critique.md): **pixel decoding and metadata extraction exist via CBFlib; miniCBF/Pilatus paths are first-class; full-CBF and packing edge cases need care.**

## Verdict (critique four-way)

Closest match: **miniCBF works (Pilatus / ADSC-PAD); other CBF variants may fail or need workarounds.**

Not “no usable CBF reader” — the tree already ships CBFlib and Fortran/C glue.

## Evidence

| Piece | Location | Role |
| ----- | -------- | ---- |
| CBFlib | [`cbf/`](../../cbf/), `cbf/lib/libcbf.a` | Decode CBF / imgCIF |
| Full/header CBF wrap | [`mosflm/cbfwrap.c`](../../mosflm/cbfwrap.c) `cbfwrap_()` | CBFlib → `cbf_int` / `cbf_int4` / `cbf_double` / `cbf_char` + `short image[]` |
| Pilatus miniCBF | [`mosflm/pilatus_c.c`](../../mosflm/pilatus_c.c), [`pilatus_pck.c`](../../mosflm/pilatus_pck.c) | Fast path; `pilatus_pck` can decode `CBF_BYTE_OFFSET` / `CBF_NONE` without full CBFlib |
| Detector keyword | [`mosflm/control.f`](../../mosflm/control.f) ~8435–8497 | `DETECTOR CBF` / Pilatus / `ADSC-PAD` set `MACHINE`, `MINICBF`, `FULLCBF` |
| Image open | `OPENODS` (amalgamation + `control.f` ~16426) | Dispatches by machine type; CBF uses CBFlib via this path |
| Shared COMMON | `CBF_PAR` in `control.f` ~472–481 | `CBF_DOUBLE(16)`, `CBF_INT(16)`, `CBF_INT4(16)`, `CBF_CHAR(16)`, `FULLCBF`, `CBFSETUP` |
| HDF5 (separate) | [`mosflm/fromhdf5.c`](../../mosflm/fromhdf5.c) | Eiger HDF5; not CBF |

## What `cbfwrap_` fills (high level)

From [`cbfwrap.c`](../../mosflm/cbfwrap.c):

- Dimensions / axis precedence → `cbf_int[]`, directions → `cbf_char[]`
- Wavelength → `cbf_double[0]`
- Detector distance → `cbf_double[1]`
- Pixel sizes (mm) → `cbf_double[2]`, `[3]`
- Gain / overload → `cbf_double[4]`, `cbf_int4[0]`
- φ start / end / range → `cbf_double[5]`, `[6]`, `[7]`
- Beam centre (mm) → `cbf_double[9]`, `[10]` (note sign flip on one axis: `-beam_centre[3]`)
- Pixels → `short image[...]` when `modeop` requests image data

`modeop == 2` is used for header-only style paths (see comments around line 526).

## Known gaps / hazards

- **Packing**: comments in `cbfwrap.c` (~601) note some full CBFs (e.g. Bruker Photon) omit packing tags.
- **Diamond CBFs**: comment ~549 about `failnez` with certain categories.
- **Hardcoded image buffer** size in the `cbfwrap_` signature (`short image[64000000]`) — large detectors / orientation mistakes are risky.
- **Multiple readers**: Pilatus path vs full `cbfwrap` can diverge in metadata conventions.
- **HDF5**: not CBF; use `fromhdf5` or external conversion (see [XDS_BRIDGE.md](XDS_BRIDGE.md)).

## Fixture for smoke tests

Symlink / copy under [`testdata/cbf/`](../testdata/cbf/README.md). Preferred sample used during this pass:

`mxflask/mxpipe/test_data/KIT-02AT05_w1_1_0001.cbf` (repo-local if present).

### Smoke result (Oct 2026)

`tools/cbf_smoke` against `testdata/cbf/sample_minicbf.cbf`:

- **Pixels OK:** `ndim1=1679 ndim2=1475`, intensity range decoded via CBFlib `cbf_get_image`.
- **Header gaps on this file:** `cbf_get_wavelength` / `cbf_construct_detector` failed (typical miniCBF layout; Pilatus uses `pilatus_c` / `pilatus_pck` for metadata). Confirms: CBFlib pixel path works; metadata routing remains detector-specific.

### `mosflm::image::CbfImageReader` (Oct 2026)

C++ reader in `diffraction-lab` links vendored `cbf/lib/libcbf.a`:

- Decodes pixels via `cbf_get_image` (CTest `cbf_image_reader`: 1475×1679, non-empty).
- Fills `ImageMetadata` from CIF when present; for miniCBF gaps, parses `#` comment
  lines inside `_array_data.header_contents` (wavelength, distance, pixel size,
  beam centre, oscillation, detector id).
- `ImageMetadata::to_detector_geometry()` maps into `mosflm::DetectorGeometry`.

**Remaining holes:** full-CBF packing edge cases; beam-centre sign conventions vs
`cbfwrap_`; no DiffractionImage (CCP4) path; no dxtbx.
