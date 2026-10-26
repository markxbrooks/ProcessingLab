# durin (adapted portions)

[durin](https://github.com/DiamondLightSource/durin) is Diamond Light
Source's HDF5 plugin for XDS. ProcessingLab does not build or link durin;
[processing-lab/src/image/Hdf5ImageReader.cpp](../../processing-lab/src/image/Hdf5ImageReader.cpp)
adapts these parts of the August 2026 source (`src/file.c`, `src/plugin.c`,
`src/filters.c`):

- NXdetector / NXdata group discovery by `NX_class`
- the single NXmx `/entry/data/data` image dataset (including virtual
  datasets) beside DECTRIS `data_NNNNNN` blocks, with blocks of unequal size
- scaling pixel size and distance by their `units` attribute
- `pixel_mask` and `detectorSpecific/pixel_mask` (any low-byte bit masks the
  pixel)
- reading bitshuffle-LZ4 chunks directly and decompressing them outside the
  HDF5 filter pipeline

Copyright (c) 2018-2019, Diamond Light Source Ltd. Distributed under the BSD
3-Clause licence in [LICENSE](LICENSE), which must accompany any
redistribution of the adapted code.
