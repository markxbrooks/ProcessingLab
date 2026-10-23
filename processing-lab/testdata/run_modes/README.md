# Input-mode fixtures (`run_modes`)

Inputs for the CTest `run_modes` (`tests/differential/test_run_modes.cpp`),
which runs `ProcessingLab` end to end on the CCP4 insulin frame
(`$CCP4/examples/data/insulin_1_001.img`, or `$MOSFLM_IMAGE`) in each input
mode. The test is skipped if the image is missing.

| File | Mode |
| ---- | ---- |
| `insulin_XDS.INP` | `run --xds-inp` (cell, symmetry, geometry; names the image) |
| `insulin.xinfo` | `run --xinfo` with no `--image` |
| `insulin_two_sweeps.xinfo` | `run --xinfo` with two sweeps: needs `--sweep`; each sweep names its own wavelength |
| `process.in` | CONTROL script mode with `--output-dir` |

`@IMAGE_DIR@` in the import files is replaced with the insulin image's
directory. The test writes the filled-in copies, and every output, under
`<build>/test_output/run_modes/`; nothing is written here or next to the
image, and the test checks that.
