# Domain C++ API vs flat C ABI

Processing kernels in `mosflm::diffraction` take **operation-specific request
objects**. The Fortran/C boundary stays flat and builds those requests.

| Operation | C++ entry | Request / result |
| --------- | --------- | ---------------- |
| Predict | `ReflectionPredictor::predict` | `PredictionRequest` → `PredictionResult` |
| Summation integrate | `SummationIntegrator::integrate` | `IntegrationRequest` |
| Profile fit | `ProfileFitter::fit` | `ProfileFitRequest` |
| Known-indexed PROCESS | `ProcessPipeline::run` | `ProcessRequest` → `ProcessResult` |
| MEAS boxes | `MeasIntegrator::integrate` | `MeasurementBoxSettings` → `MeasObservation` |
| Standard profile | `StandardProfileFitter::fit` | `StandardProfileSettings` → `ProfileObservation` |
| Post-refinement | `PostRefiner::refine` | observations + settings → `PostRefinementResult` |
| CONTROL script | `ControlSession::run_script` | `.in` keywords → `ControlGoResult` / PROCESS |
| LP factor | `lp_factor(detector, x, y)` | small helper (no request object) |

## Prediction

- `ScanRange` — physical φ start/end (still when equal).
- `PredictionSettings` — `dmin_A`, Ewald tolerance, **algorithm** `phi_step_deg`, mosaic + `divh_deg`/`divv_deg` (REEK E0 gate widening)
  (search increment, not oscillation width).
- `PredictionRequest` — crystal, orientation, detector, scan, settings.
- `PredictionResult::reflections` — predicted spots.

## Integration

`IntegrationRequest` / `ProfileFitRequest` hold **references** to the image,
detector, and a `std::span` of predictions plus settings. Referenced objects
must outlive the call (no giant owned context blob).

## PROCESS (Phase D / D1–D5)

`ProcessRequest` / `ProcessResult` own the known-indexed pipeline without
Fortran GENERATE / REEK / MEAS / INTEG3 / PROCESS:

1. Optional `PostRefiner` when `enable_postrefinement` (default **false**)
2. `ReflectionPredictor` → in-memory spot table (`predictions`)
3. `StandardProfileFitter` (default) — MEAS boxes + learned Rossmann standard
   profile (`I = ΣP·(OD−bg)·ΣP / ΣP²`); optional `MeasIntegrator` /
   `SummationIntegrator`
4. `MtzAdapter::write_simple` via `ProcessPipeline::write_mtz`

`ProcessPipeline::known_indexed_insulin` uses geometry from
`testdata/images/known_indexed/geometry.txt` / `ipmosflm.in` (SEPARATION
0.57 mm). CTests: `insulin_cpp_process`, `insulin_isi`, `insulin_profile`.

## CONTROL (Phase D / D5)

`ControlSession` parses a MOSFLM-like keyword subset (CELL, SYMMETRY, BEAM,
DISTANCE, WAVELENGTH, PIXEL, MATRIX/AMAT, PROCESS, GO, EXIT, …) and on GO
builds a `ProcessRequest` + runs `ProcessPipeline`. Script:
`testdata/images/known_indexed/cpp_process.in`. CLI:
`mosflm_cpp_cli run <script.in> [--image img] [mtz]`. Thin default binary:
`ipmosflm` (Phase D / D6; no amalgamation). CTests: `control_session`,
`insulin_control`, `ipmosflm_thin`. Not a line port of Fortran CONTROL
(no AUTOINDEX / STRATEGY / iMosflm sockets).

## Post-refinement (Phase D / D4)

`PostRefiner` solves spindle-+X Ewald φ, then coordinate-descent on missetting
(`dφx/y/z`), isotropic `cell_scale`, and mosaic from MAD of φ residuals.
Optional ``refine_against_detector`` adds XDET/YDET positional cost via
``predict_xdet_ydet_px`` (POSTREFL parameter set; not full IDXREF / multi-image
NADD). CTests: `post_refinement`, `insulin_postref`. Keep off by default in
PROCESS until caller workflows opt in.
Multi-region ``NMASKS``: ``StandardProfileSettings::n_bins_x/y`` (PROCESS default
``5×5``); sparse regions fall back to the global learned profile.

## Detector geometry

``DetectorGeometry::lab_to_detector_mm`` projects onto ``z = distance``, then
applies MOSFLM ``MMTOPX`` tilt/twist: default old ``DXTOF = XTOFRA + TILT·X +
TWIST·Y``, or ``use_rotmat2`` (NEWT) ray–plane scale from ``rotmat2_detnor``.
TILT/TWIST fields are printed milliradian integers (``ITILT``/``ITWIST``).

## Flat ABI

[`C_ABI.md`](C_ABI.md) and `MosflmCApi.h` remain the interoperability surface.
Adapters in `MosflmCApiValidation.cpp` / `LegacyMosflmAdapter` assemble domain
requests from POD arguments.
