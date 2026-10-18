#pragma once

#include "mosflm/image/DiffractionImage.hpp"

#include <filesystem>

namespace mosflm::image {

/// Read a diffraction image by extension (``.img``/``.smv`` → SMV, ``.cbf`` → CBF).
///
/// HDF5 / NeXus (``.h5`` / ``.nxs``) are recognized but not yet decoded —
/// throws with a clear Phase E message.
[[nodiscard]] DiffractionImage read_diffraction_image(
    const std::filesystem::path& path
);

}  // namespace mosflm::image
