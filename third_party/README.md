# Third-party code

Code from other projects, kept apart from ProcessingLab's and MOSFLM's own
sources. Each folder carries its upstream licence; keep it with the code.
Note any local changes in this table so upgrades do not lose them.

| Folder | Upstream | Licence | Used by |
| ------ | -------- | ------- | ------- |
| [`bitshuffle/`](bitshuffle/) | Bitshuffle, Kiyoshi Masui | MIT ([LICENSE](bitshuffle/LICENSE)) | HDF5 reading in ProcessingLab (CMake) and ipmosflm (`mosflm/Makefile`) |
| [`lz4/`](lz4/) | LZ4, Yann Collet; `h5zlz4.c` from the DECTRIS HDF5 plugin, modified for static linking by Takanori Nakane | BSD 2-Clause ([LICENSE](lz4/LICENSE)) | As bitshuffle |
| [`cbf/`](cbf/) | CBFlib, Paul Ellis and Herbert J. Bernstein | GPL or LGPL ([COPYING](cbf/COPYING)) | CBF images in ProcessingLab and ipmosflm (`libcbf.a`, built by `make cbf-lib`) |
| [`jpg/`](jpg/) | Independent JPEG Group libjpeg, release 6b | IJG ([README](jpg/README)) | Nothing at present: `configure` links the system libjpeg |
| [`durin/`](durin/) | durin, Diamond Light Source | BSD 3-Clause ([LICENSE](durin/LICENSE)) | Licence for code adapted into `processing-lab/src/image/Hdf5ImageReader.cpp`; see [durin/README.md](durin/README.md) |

Local changes: none recorded yet.
