#pragma once

#include "mosflm/image/DiffractionImage.hpp"

#include <filesystem>

namespace mosflm::image {

/// Read a diffraction image by extension.
///
/// - ``.img`` / ``.smv`` → SMV
/// - ``.cbf`` → CBF
/// - ``.h5`` / ``.nxs`` / ``.hdf5`` → HDF5/NeXus (Eiger master / DLS NXmx)
///
/// ``frame`` is 1-based and used for multi-image HDF5 masters (CONTROL
/// ``PROCESS`` first image). SMV/CBF ignore ``frame``.
[[nodiscard]] DiffractionImage read_diffraction_image(
    const std::filesystem::path& path,
    int frame = 1
);

}  // namespace mosflm::image
