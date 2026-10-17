# First Fortran routine to replace: `PREDICT`

## Why this routine

The vertical-slice plan is: **known orientation → predict spots → summation integrate**.
`PREDICT` is the first diffraction-science kernel after geometry that produces
observable outputs (predicted detector coordinates / φ ranges) without requiring
a full command interpreter port.

Autoindexing (DPS) is deferred; REFIX is unavailable in the public tree.

## Location

| Item | Value |
| ---- | ----- |
| Primary | `mosflm/mosflm_all_ip_inc_noX.f` ≈ line **171805** (`SUBROUTINE PREDICT(PSI,IFLAGPR)`) |
| X11 twin | `mosflm/mosflm_all_ip_inc.f` ≈ line **180578** |
| Closely related | `DSTAR` (`unoptimizable.f:17692`), `ROTMAT` (`unoptimizable.f:37480`), `CLCALC` / `MATSET` |

Prefer the **no-X** amalgamation for headless differential tests.

## Callers (high level)

`PREDICT` is driven from the integration / measurement path (`PROCESS` /
`MEAS` / control flow after CELL + MATRIX are set). Exact call graph is dense
inside the amalgamation; when porting, extract:

1. Inputs read from COMMON / arguments at entry.
2. Side effects written back (predicted spot lists, φ beg/end, detector x/y).
3. Downstream consumers (`MEAS`, display overlays — ignore display initially).

## I/O contract (to freeze before coding)

**In (minimum):**

- Unit cell + orientation (UB / MOSFLM MATRIX).
- Wavelength, crystal-to-detector distance, beam centre, pixel size.
- Oscillation start/end (or single φ for stills).
- Resolution / packing limits used by MOSFLM for this run.

**Out (minimum for regression):**

- List of predicted reflections: `h k l`, detector `x y` (mm or pixels),
  φ begin/end (or still flag), full/partial classification.
- Count of predictions in the image window.

Document the precise MOSFLM arrays once a single scripted `ipmosflm` job is
captured under `testdata/oracle/`.

## Regression criteria

| Check | Tolerance (starting point) |
| ----- | -------------------------- |
| Same (h,k,l) set on a still/image | Exact set match after sorting |
| Detector coordinates (mm) | ≤ 0.05 mm or ≤ 0.5 pixel |
| φ mid / width for partials | ≤ 0.05° |
| Spot count on detector | Exact or ±1 for edge cases |

Do **not** require bit-identical floats.

## Prerequisites before swapping callers

1. `mosflm::UnitCell` / `OrientationMatrix` / `DetectorGeometry` unit tests green.
2. Oracle `bin/ipmosflm` job on `testdata/` producing a frozen prediction dump.
3. C++ `ReflectionPredictor` API matching the I/O contract above.
4. Differential test comparing C++ vs oracle dump.
5. Only then route one caller through `LegacyMosflmAdapter` → C++.

## Non-goals for the first replacement

- Profile fitting (`INTEG3` / `INTEG4`).
- Post-refinement.
- iMosflm socket / X11 overlays.
- Command-language compatibility (`control.f`).
