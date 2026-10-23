I agree. For the first version, ProcessingLab should have a minimal, image-first CLI. Let the image headers supply experiment metadata, with optional `.xinfo` or `XDS.INP` files to override or supplement it. Avoid requiring a YAML configuration file until there is a demonstrated need.

I'd start with this interface.

## 1. Proposed CLI

```
# Simplest run: discover metadata and process the dataset
ProcessingLab run --image /data/sweep_master.h5

# Specify the output directory explicitly
ProcessingLab run --image /data/sweep_master.h5 \
    --output-dir results/lysozyme

# Limit processing to an image range
ProcessingLab run --image /data/sweep_master.h5 \
    --range 1:180

# Use XDS.INP to override or supplement image metadata
ProcessingLab run --image /data/sweep_master.h5 \
    --xds-inp /data/XDS.INP

# Use Xia2 .xinfo instead
ProcessingLab run --image /data/sweep_master.h5 \
    --xinfo /data/automatic.xinfo

# Combine both sources, with explicit precedence
ProcessingLab run --image /data/sweep_master.h5 \
    --xds-inp /data/XDS.INP \
    --xinfo /data/automatic.xinfo \
    --range 1:180 \
    --output-dir results/lysozyme
```

I'd support multiple `--image` arguments from the start, even if the first implementation processes one sweep at a time:

```
ProcessingLab run \
    --image /data/sweep1_master.h5 \
    --image /data/sweep2_master.h5
```

For CBF or SMV sequences, accept an individual image or a filename template, depending on the reader's capabilities. HDF5 masters should be treated as a single image source, not as one file per frame.

## 2. Output directory: timestamped by default

I would make automatic output-directory creation the default. Never silently reuse an existing output directory.

Default: unique timestamped directory

```
ProcessingLab run --image sweep_master.h5
```

Example output:

```
processing_20261022_211530/
```

Explicit directory

```
ProcessingLab run --image sweep_master.h5 \
    --output-dir results/lysozyme
```

If the directory already exists, fail clearly unless a safe resume or overwrite policy is explicitly selected.

Optional named run

```
ProcessingLab run --image sweep_master.h5 \
    --output-dir results/lysozyme \
    --run-name trial_02
```

A useful addition later, but not essential for the first pass.

I'd use a timestamp with second resolution, plus a collision check. If two jobs start simultaneously, the second must not overwrite the first; it should generate a unique suffix or report the conflict.

## 3. Metadata precedence

This is the most important behavior to define correctly. I suggest the following precedence, from lowest to highest priority:

| Priority | Source               | Role                                                     |
| -------- | -------------------- | -------------------------------------------------------- |
| 1        | Image/HDF5 metadata  | Default geometry, wavelength, detector and scan metadata |
| 2        | `XDS.INP`            | Explicit geometry, cell, symmetry and scan settings      |
| 3        | `.xinfo`             | Sweep selection, image paths, frame range and wavelength |
| 4        | Explicit CLI options | Image range and output location                          |

The `.xinfo` and `XDS.INP` formats overlap, so precedence should apply per field, not by replacing one entire configuration with another. In particular, `.xinfo` should not erase detector geometry obtained from the image headers merely because it lacks those fields.

I would also make the final resolved configuration visible in the log:

```
Experiment metadata
  Image source:  /data/sweep_master.h5
  Wavelength:    0.9795 Å  [HDF5 header]
  Distance:      200.0 mm  [HDF5 header]
  Cell:          57.98 58.04 148.16 ... [XDS.INP]
  Frame range:   1–180    [CLI]
  Output:        processing_20261022_211530/
```

If the sources disagree, record which value won and warn about scientifically meaningful conflicts. Do not silently reconcile contradictory geometry.

One nuance: a `.xinfo` file can describe multiple sweeps. For the first release, either require an unambiguous sweep or provide a simple `--sweep` selector rather than silently choosing the first sweep.

## 4. Keep the initial option set small

I'd implement only these options initially:

| Option               | Required?                              | Purpose                                          |
| -------------------- | -------------------------------------- | ------------------------------------------------ |
| `--image PATH`       | Yes, unless supplied by an import file | Image or HDF5 master; repeatable                 |
| `--range FIRST:LAST` | No                                     | Inclusive frame range                            |
| `--output-dir DIR`   | No                                     | Explicit output directory; otherwise timestamped |
| `--xds-inp PATH`     | No                                     | Import XDS configuration                         |
| `--xinfo PATH`       | No                                     | Import Xia2 experiment description               |
| `--help`             | No                                     | Display usage                                    |
| `--version`          | No                                     | Display build version                            |

I would omit `--cell`, `--space-group`, `--beam`, `--distance`, `--wavelength`, `--threads` and algorithm-specific options from the initial public CLI. Those can remain available through CONTROL until a genuine use case warrants exposing them.

For a first implementation, I would also avoid a separate `run` subcommand if you prefer to keep the executable invocation even simpler. Both of these are reasonable:

```
ProcessingLab --image sweep_master.h5 --range 1:180
```

```
ProcessingLab run --image sweep_master.h5 --range 1:180
```

My preference is the second, because it leaves room for `inspect`, `index` and `integrate` without changing the interface later.

## 5. One implementation detail worth getting right

Your current import logic already has a precedence order for `--xds-inp`, `--xinfo` and `--xparm`. I'd preserve that behavior where it is appropriate, but distinguish imported values from values read directly from image headers. The new CLI should resolve the complete experiment model once, validate it, and then pass it to the processing engine.

That gives you a small first release with a clear contract:

- Specify images.
- Optionally restrict the frame range.
- Optionally import XDS or Xia2 metadata.
- Automatically create a fresh output directory.
- Resolve metadata transparently and report its provenance.
- Execute the existing processing pipeline.

That is enough to establish a clean CLI without prematurely designing a large configuration language.