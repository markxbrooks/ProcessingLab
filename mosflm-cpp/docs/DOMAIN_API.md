# Domain C++ API vs flat C ABI

Processing kernels in `mosflm::diffraction` take **operation-specific request
objects**. The Fortran/C boundary stays flat and builds those requests.

| Operation | C++ entry | Request / result |
| --------- | --------- | ---------------- |
| Predict | `ReflectionPredictor::predict` | `PredictionRequest` → `PredictionResult` |
| Summation integrate | `SummationIntegrator::integrate` | `IntegrationRequest` |
| Profile fit | `ProfileFitter::fit` | `ProfileFitRequest` |
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

## Flat ABI

[`C_ABI.md`](C_ABI.md) and `MosflmCApi.h` remain the interoperability surface.
Adapters in `MosflmCApiValidation.cpp` / `LegacyMosflmAdapter` assemble domain
requests from POD arguments.
