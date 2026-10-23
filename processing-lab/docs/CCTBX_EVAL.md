# CCTBX selective evaluation (Bravais / change-of-basis)

Status: **defer hard dependency**; **adopt as optional reference validator**.

## Context

Critique asked whether CCTBX (`sgtbx` / `uctbx`) should replace or sit beside
Gemmi for crystallographic infrastructure. Local checkout:
[`/Users/brooks/projects/cctbx_project`](../../../cctbx_project) (sibling of
`ipmosflm`). **dxtbx is not present** there and remains out of scope; CBF stays
on vendored CBFlib.

## Build / run cost

| Approach | Effort on this machine | Notes |
| -------- | ---------------------- | ----- |
| C++ `libtbx` link into `processinglab_geometry` | **High** | No ready `cctbx_build` found under `~/projects`; full bootstrap is multi-hour and would pull Boost/scitbx into every target |
| Phenix `cctbx.python` | **Low** | `/Applications/phenix-1.20.1-4487/build/bin/cctbx.python` (Python **2.7**) already imports `cctbx.uctbx` / `cctbx.sgtbx`; POC script is Py2/Py3-safe |
| Conda `cctbx-base` | Medium | Documented alternative; not required for the POC |

**Decision:** first-pass evaluator is a **standalone Python POC**
([`tools/cctbx_cell_poc.py`](../tools/cctbx_cell_poc.py)). CMake keeps
`MOSFLM_WITH_CCTBX=OFF` by default and does **not** link CCTBX into
`processinglab_geometry`.

```bash
# From processing-lab/
/Applications/phenix-1.20.1-4487/build/bin/cctbx.python tools/cctbx_cell_poc.py
```

Expect `CCTBX_POC_OK`.

## Numerical agreement (POC, Oct 2026)

Cell `(58, 72, 90, 90, 90, 90)`:

| Quantity | CCTBX | `mosflm::UnitCell` / INVCEL | Gemmi (via `test_gemmi_crystal`) |
| -------- | ----- | --------------------------- | -------------------------------- |
| Volume | matches within `1e-4` | reference | matches MOSFLM within `1e-4` |
| Reciprocal lengths / angles | matches within `1e-8` / `1e-4` | reference | `a*`, `c*` match INVCEL |

Change-of-basis:

- Identity `a,b,c` leaves parameters unchanged.
- DPS lattice-character matrix #2
  (`{1,-1,0; -1,0,1; -1,-1,-1}` from `bravais_transform.c`) applied via
  `sgtbx.change_of_basis_op` scales volume by `|det|` as expected.

## API surface used

- `cctbx.uctbx.unit_cell` — parameters, volume, reciprocal, `change_basis`
- `cctbx.sgtbx.change_of_basis_op` / `rt_mx` — integer CoB from DPS `tm[]`
- `scitbx.matrix.sqr` — determinant / matrix construction

Not used (and not needed for Bravais validation): `dxtbx`, `rstbx`, full
`cb_op` databases beyond explicit matrices.

## Decisions per subsystem

| Subsystem | Decision | Rationale |
| --------- | -------- | --------- |
| Default domain layer (`UnitCell`, MTZ, space groups) | **Keep Gemmi** | Already in-tree; CTest `gemmi_crystal` green; zero CCTBX link cost |
| DPS Bravais algorithm | **Keep port / façade** (`mosflm::indexing::BravaisClassifier`) | Proven MOSFLM path; critique: wrap, do not rewrite yet |
| CCTBX as Bravais / CoB validator | **Adopt (optional)** | Independent metric check when Phenix/`cctbx.python` available |
| CCTBX as hard dependency of `processing-lab` | **Defer / reject for day one** | `libtbx` bootstrap cost; no dxtbx locally; Gemmi covers infra |
| dxtbx image models | **Out of scope** | Missing from local `cctbx_project`; CBF → CBFlib |

## CMake hook

```cmake
option(MOSFLM_WITH_CCTBX "Optional CCTBX validator targets" OFF)
```

When `ON` in a future pass, add a non-default target that shells to
`cctbx_cell_poc.py` or links a tiny C++ harness — never into the default
`processinglab_geometry` library.

## Next science kernel

After Bravais façade + `CbfImageReader`:
[`FIRST_ROUTINE.md`](FIRST_ROUTINE.md) (`PREDICT` / `ReflectionPredictor`).
