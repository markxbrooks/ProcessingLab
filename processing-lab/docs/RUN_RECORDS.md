# Run records: `input.json` and `process.json` (schema version 2)

Every `ProcessingLab run`, and every CONTROL script run with `--output-dir`,
writes two JSON records into its output directory:

- **`input.json`:** what the run was asked to process, resolved before
  processing starts. It covers images, imports, every experiment value with
  its source, and the effective settings. It is written once.
- **`process.json`:** what happened. It covers stage status, the indexing
  solution, the refined geometry, integration counts, quality statistics,
  warnings, errors and outputs. It is rewritten atomically after every
  change.

`process.md` is rendered from the two records and is never parsed.

## Compatibility

- **Versioning:** both records carry `schema_version`. A change that renames
  or removes a field, or changes its meaning or units, increments the
  version.
- **Unknown fields:** readers must ignore them, so new fields can be added
  without a version change.
- **Older runs:** `run --from` reads version 1 and version 2 runs. It refuses
  a newer version, rather than guessing at its layout.
- **Changes in version 2:**

  | Version 1 | Version 2 |
  | --------- | --------- |
  | `started` / `finished` / `created`, `-0700` offsets | `started_at` / `finished_at` / `created_at`, RFC 3339 (`-07:00`), plus `duration_s` |
  | `"unit": "angstrom/deg"` | `length_unit: "angstrom"`, `angle_unit: "degree"` |
  | `space_group_number` | `space_group: {number, symbol}` |
  | `beam_centre_mm: [x, y]` | `beam_centre_mm: {x, y}` |
  | `tilt_twist_mrad: [tilt, twist]` | `detector_tilt_twist_mrad: {tilt, twist}` |
  | `refinement` (final values, `scope` free text) | `experiment` (the final model, with conventions); `refinement` (`stage`, what was refined or fixed) |
  | Absolute output paths | Paths relative to `run_directory` |
  | `integration.predictions` / `reflections` / `profile_fitted` / … | `integration.counts`, defined below |
  | `quality.note`, three counters | `scaling_status`, `outlier_rejection_status`, range, shell definition, `quality.counts` |
  | Stage keys present only when set | Every stage has every key; unset ones are `null` |
  | — | `validation` |

## Conventions

- **Units:** these are in the field name (`_mm`, `_deg`, `_mrad`,
  `_angstrom`, `_s`) or in an explicit `*_unit` field. Unit cells give
  lengths in Å and angles in degrees.
- **Timestamps:** RFC 3339 local time with its UTC offset. `duration_s` is
  measured on a monotonic clock, so it doesn't depend on the timestamps.
- **Paths:** output paths are relative to `run_directory`, the absolute
  directory where the run was written. Image and import paths in
  `input.json` are recorded as given, made absolute, and never rewritten.
- **Missing values:** a missing value is `null`, never a sentinel number.
- **Orientation matrix (`mosflm_amat`):** MOSFLM [A] = λ·U·B.
  - **Rows:** the MOSFLM lab axes: X along the beam, Y, and Z the spindle.
  - **Columns:** a\*, b\*, c\*.
  - **Scale:** `orientation.wavelength_angstrom` (λ).
  - **`invertx`:** the detector handedness the matrix was derived in, as the
    CONTROL `INVERTX` keyword sets it. ProcessingLab's own indexing and XPARM
    files use `false`; a typed MOSFLM matrix uses `true`.
  - **`orientation.unit_cell`:** the cell the matrix implies. In a final model
    it equals `crystal.unit_cell`.
- **Beam centre:** `{x, y}` in mm along ProcessingLab's detector axes, as
  header pixel coordinates × pixel size.
  - **CBF and HDF5:** x is the image's fast axis and y its slow axis.
  - **ADSC SMV:** the image is stored transposed, as MOSFLM reads it, and the
    header `BEAM_CENTER_X/Y` is used as it stands.
  - **Same frame:** the CONTROL `BEAM` keyword and `XDS.INP` `ORGX×QX`,
    `ORGY×QY`.
- **Which matrix to carry forward:** `experiment.orientation`, the final one.
  `indexing.selected_solution.orientation_matrix_mosflm` is the autoindexing
  matrix, kept as history; it implies the unconstrained indexing cell.

## `process.json`

### Run

| Field | Meaning |
| ----- | ------- |
| `status` | `running` → `completed` or `failed`. A record left at `running` means the process died. |
| `started_at`, `finished_at`, `duration_s` | Run timing; the last two are `null` while running. |
| `run_directory` | Absolute path of the run directory. |
| `validation` | `{status, checks}`. Always `not_assessed` with no checks for now. `completed` means the program ran to the end, not that the solution or data are correct. |
| `error` | The error that failed the run, or `null`. |
| `warnings` | Run-level warnings, starting with those from metadata resolution. |
| `outputs` | `[{kind, path}]` for files that exist. `path` is relative to `run_directory`. |

### Stages

`stages` lists the planned steps in order. Each is an object with the same
keys:

| Key | Meaning |
| --- | ------- |
| `name` | `findspots`, `autoindex`, `process` or `strategy` |
| `status` | `pending` → `running` → `completed` or `failed`; `skipped` if the run finished without reaching it |
| `started_at`, `finished_at`, `duration_s` | Stage timing, or `null` |
| `message` | Human-readable summary, or `null`. Facts are in `results`; don't parse this. |
| `results` | Machine-readable summary, or `null`. `findspots`: `{spots}`. `autoindex`: `{lattice, indexed_fraction, rmsd_hkl}`. `process`: `{integrated_reflections}`. |
| `error` | The stage's error, or `null` |

A completed stage means it ran. It doesn't mean its result is right: a
completed `autoindex` can still have chosen the wrong lattice.

### `indexing.selected_solution`

The autoindexing result, kept as history. Its `unit_cell` is the
unconstrained indexing estimate, and is not updated by refinement.

| Field | Meaning |
| ----- | ------- |
| `lattice`, `lattice_character` | Bravais lattice (e.g. `tP`) and DPS lattice character |
| `space_group` | `{number, symbol}`, the symmetry assumed for refinement |
| `unit_cell` | `{a, b, c, alpha, beta, gamma, length_unit, angle_unit}` |
| `orientation_matrix_mosflm` | 3×3 rows of the autoindexing matrix, in ProcessingLab's internal MOSFLM-style convention (not yet specified as an interchange format) |
| `spots_used`, `indexed_fraction`, `rmsd_hkl` | Spots used, the fraction indexed and the RMS deviation from integer hkl |

`indexing.solutions` lists every lattice solution AUTOINDEX considered (the
penalty table in `process.log`), best first:

| Field | Meaning |
| ----- | ------- |
| `rank` | 1 = lowest penalty |
| `character` | Lattice character 1–44; `AUTOINDEX … CHARACTER n` selects this solution |
| `lattice` | Bravais lattice, e.g. `tP` |
| `penalty` | DPS penalty (lower is better) |
| `unit_cell` | Conventional cell for that lattice |
| `metric_deviation_pct` | Departure from the lattice's symmetry-constrained metric (%), or `null` |
| `space_groups` | Possible space groups, lowest symmetry first |
| `selected` | Whether AUTOINDEX chose it |

`indexing.lastgood_rank` is the rank of the last solution before the largest
penalty jump (MOSFLM's LASTGOOD).

### `experiment`

The run's final model, after the last periodic refinement of `PROCESS`. It
is `null` until `PROCESS` has run. Every value is stated once, here, and is
what `run --from` carries forward. Values that CONTROL didn't hold (for
example a wavelength from the image header) come from `input.json`, never
from built-in defaults.

| Field | Meaning |
| ----- | ------- |
| `crystal.unit_cell` | Final cell (`length_unit` / `angle_unit`), constrained to the space group's metric |
| `crystal.space_group` | `{number, symbol}` |
| `wavelength_angstrom` | Wavelength |
| `detector.distance_mm` | Distance |
| `detector.beam_centre_mm` | `{x, y}` (see Conventions) |
| `detector.pixel_size_mm` | `{x, y}` |
| `detector.tilt_twist_mrad` | `{tilt, twist}`; `GXPARM.XDS` doesn't carry these |
| `detector.size_px` | `{fast, slow}` from the header |
| `scan` | `{frames: {first, last}, oscillation_deg}` |
| `mosaicity_deg` | Final mosaicity |
| `orientation` | `{convention: "mosflm_amat", matrix, wavelength_angstrom, invertx, unit_cell}` (see Conventions), or `null` |
| `conventions` | The two conventions above, in words, so the record stands alone |

### `refinement`

What `PROCESS` did to each parameter: the policy it ran with.

| Field | Meaning |
| ----- | ------- |
| `stage` | `process` |
| `parameters` | `orientation`, `unit_cell`, `detector_distance`, `beam_centre`, `detector_tilt_twist`: each `refined` or `fixed`. `mosaicity`: `estimated` from the data or `fixed` |
| `unit_cell_constraint` | `space_group_metric` when the cell is held to the space group's metric, otherwise `none` |

The policy the run was asked for is `experiment.refinement_policy` in
`input.json`.

## Second-round processing (`run --from`)

`run --from` builds the next run's starting model from the previous run's
`process.json` `experiment` block. It does that only after checking the
block. It refuses:

- a run that isn't `completed`, or never ran `PROCESS`;
- a model missing any geometry field;
- a cell that its space group's crystal system forbids;
- an orientation convention other than `mosflm_amat`;
- a matrix that implies a cell other than `crystal.unit_cell`. The check
  recomputes that cell rather than trusting the stored copy, and allows
  0.2% on lengths and 0.2° on angles;
- a wavelength or pixel size that disagrees with that run's `input.json`.

The previous run is never modified. The new run's `input.json` records it
under `inputs.imports.previous_run`, and tags every inherited value
`previous_run`. Values the command line overrides (`--range`, `--cell`,
`--spacegroup`) keep both, under `overrides`.

The refinement policy is explicit, recorded in `input.json`
(`experiment.refinement_policy`) and in `process.json`
(`refinement.parameters`). `--refine` sets it; parameters it doesn't name
are held fixed.

| `--refine` | Effect |
| ---------- | ------ |
| not given, with `--from` / `--xparm` | `orientation,mosaicity`: refine the orientation and re-estimate the mosaicity; cell, distance, beam and tilt/twist stay at their previous values |
| not given, otherwise | Everything is refined |
| `orientation,cell,distance,beam,tilt_twist,mosaicity` (any subset), `all`, `none` | Exactly those. Leaving out `mosaicity` holds it at the previous value (`MOSAIC`), which needs `--from` |

Runs without an `experiment` block (version 1, or runs that stopped before
`PROCESS`) carry the orientation through `GXPARM.XDS`, and the tilt/twist
through `refinement.tilt_twist_mrad`, as before.

### `integration`

| Field | Meaning |
| ----- | ------- |
| `frames` | `{first, last}` integrated |
| `peak_box_px` | Measurement box width (px) |
| `counts` | See below |
| `defpix` | Pixels DEFPIX classified before integrating, or `null` with `DEFPIX OFF`: `frames_sampled`, `block_px` (shadow-test block edge; `null` with fewer than 3 frames), `defects_checked` (false when the frames span under 5°, too little to tell a defective or hot pixel from a spot) and `pixels` |

`defpix.pixels` counts every pixel once: `trusted`, `detector_masked` (gaps
and the detector's bad-pixel mask), `untrusted_region`, `excluded_resolution`,
`shaded`, `defective` and `hot`. Integration ignores all but `trusted`.
Its control image is listed in `outputs` with kind `defpix_image`; the
settings are in `input.json` `settings.integration.defpix`.

The `counts` come from two populations, which shouldn't be compared
directly:

| Count | Population |
| ----- | ---------- |
| `predicted_observations` | Predicted observations, summed over every frame, before partials are summed. A reflection spread over several frames counts once per frame. |
| `profile_fitted_observations` | Those predicted observations whose profile fit succeeded |
| `overlapped_observations` | Those predicted observations flagged as overlapped |
| `omitted` | Predicted observations not integrated, plus summed reflections dropped because they lie outside the φ range with less than half their rocking curve recorded |
| `integrated_reflections` | Reflections written to the MTZ, after partials are summed |

`profile_fitted_observations` can therefore exceed `integrated_reflections`.

### `quality`

Descriptive statistics of the integrated reflections. They are **unscaled**
(`scaling_status`) and **not outlier-rejected**
(`outlier_rejection_status`), so they show how integration went and don't
decide anything. In particular, no resolution cutoff is applied from them.

| Field | Meaning |
| ----- | ------- |
| `space_group` | `{number, symbol}` used to merge symmetry equivalents and count possible reflections |
| `resolution_range_angstrom` | `{low, high}` actually used |
| `min_i_over_sigma` | Observations below this I/σ are left out |
| `shell_definition` | `{method: "equal_width_in_inverse_d", count, d_min_angstrom}`: shells of equal width in 1/d; each shell's `d_min_angstrom` is its high-resolution edge |
| `shells`, `total` | Per-shell and overall `observed`, `unique`, `possible`, `completeness_pct`, `r_obs_pct`, `r_meas_pct`, `mean_i_over_sigma`, `cc_half_pct`, `compared` (`null` where not defined) |
| `counts` | See below |

`counts` partition the integrated reflections. Each row is assigned to the
first filter it fails, so the six filter counts add up to
`integrated_reflections`:

| Count | Rows |
| ----- | ---- |
| `integrated_reflections` | All MTZ rows (the same number as `integration.counts.integrated_reflections`) |
| `rejected_misfit_observations` | Flagged by integration (any flag except the negative-profile bit) |
| `invalid_sigma_observations` | σ ≤ 0, or I or σ not finite |
| `systematic_absence_observations` | Systematically absent in `space_group` |
| `below_min_i_over_sigma_observations` | I/σ below `min_i_over_sigma` |
| `outside_resolution_observations` | d outside `resolution_range_angstrom` |
| `accepted_observations` | The rest: the observations the statistics use |
| `unique_reflections` | Distinct asymmetric-unit hkl among the accepted observations |

The statistics are computed on the accepted observations, grouped by
asymmetric-unit hkl under `space_group`. A different space group therefore
changes them, even for the same images. For a unique reflection with n
observations, ⟨I⟩ is their 1/σ²-weighted mean:

- **`r_obs_pct`:** 100 · Σ|Iᵢ − ⟨I⟩| / Σ|Iᵢ|, over reflections with n ≥ 2.
  `compared` counts the observations included.
- **`r_meas_pct`:** the same, with each reflection's deviations scaled by
  √(n/(n−1)).
- **`mean_i_over_sigma`:** the mean, over unique reflections, of the merged
  ⟨I⟩/σ(⟨I⟩), where σ(⟨I⟩) = (Σ 1/σᵢ²)^(−1/2).
- **`cc_half_pct`:** 100 × the Pearson correlation between the means of
  random half-datasets, over reflections with n ≥ 2. It is `null` when
  nothing was compared.

## `input.json`

| Field | Meaning |
| ----- | ------- |
| `created_at` | RFC 3339 timestamp |
| `mode`, `run_mode`, `parameter_source` | `run` or `script`; `RunMode.*`; where the parameters came from, in words |
| `command_line` | The arguments, as given |
| `inputs.images` | `[{path, source, format, template?, frames_available}]` |
| `inputs.imports` | `xds_inp`, `xinfo`, `xinfo_sweep`, `xparm`, `previous_run`, each a path or `null` |
| `experiment` | Header facts `detector_size_px`, `detector_id`, `first_frame` (first frame on disk) and `phi_start_deg` (its φ), the `refinement_policy`, and each value as `{value, unit (or length_unit / angle_unit), source, overrides: [{value, source}]}`; `source` is `image_header`, `xds_inp`, `xinfo`, `xparm`, `previous_run`, `control` or `command_line` |
| `settings` | Effective FINDSPOTS / AUTOINDEX / integration / refinement settings, including the parameters refined |
| `warnings` | Disagreements between sources (the higher-priority value is used) |
| `control` | The CONTROL text the run executed |
