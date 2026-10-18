# Domain C++ API vs flat C ABI

Processing kernels in `mosflm::diffraction` take **operation-specific request
objects**. The Fortran/C boundary stays flat and builds those requests.

| Operation | C++ entry | Request / result |
| --------- | --------- | ---------------- |
| Predict | `ReflectionPredictor::predict` | `PredictionRequest` → `PredictionResult` |
| Summation integrate | `SummationIntegrator::integrate` | `IntegrationRequest` |
| Profile fit | `ProfileFitter::fit` | `ProfileFitRequest` |
| Known-indexed PROCESS | `ProcessPipeline::run` | `ProcessRequest` → `ProcessResult` |
| LP factor | `lp_factor(detector, x, y)` | small helper (no request object) |

## Prediction

- `ScanRange` — physical φ start/end (still when equal).
- `PredictionSettings` — `dmin_A`, Ewald tolerance, **algorithm** `phi_step_deg`
  (search increment, not oscillation width).
- `PredictionRequest` — crystal, orientation, detector, scan, settings.
- `PredictionResult::reflections` — predicted spots.

## Integration

`IntegrationRequest` / `ProfileFitRequest` hold **references** to the image,
detector, and a `std::span` of predictions plus settings. Referenced objects
must outlive the call (no giant owned context blob).

## PROCESS (Phase D / D1)

`ProcessRequest` / `ProcessResult` own the known-indexed pipeline without
Fortran GENERATE / REEK / PROCESS:

1. `ReflectionPredictor` → in-memory spot table (`predictions`)
2. `SummationIntegrator` (optional MOSFLM detector axes for SMV) → `integrated`
3. `MtzAdapter::write_simple` via `ProcessPipeline::write_mtz`

`ProcessPipeline::known_indexed_insulin` loads geometry / orientation from
`testdata/images/known_indexed/geometry.txt` conventions onto an SMV
`DiffractionImage` (see `SmvImageReader`). CTest `insulin_cpp_process` gates
ASU-aware HKL recall ≥ 0.85 vs the Fortran oracle MTZ (no `mosflm_predict_hook`).

## Flat ABI

[`C_ABI.md`](C_ABI.md) and `MosflmCApi.h` remain the interoperability surface.
Adapters in `MosflmCApiValidation.cpp` / `LegacyMosflmAdapter` assemble domain
requests from POD arguments.
