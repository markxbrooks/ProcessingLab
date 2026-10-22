#pragma once

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/image/DiffractionImage.hpp"

#include <filesystem>
#include <optional>

namespace mosflm::indexing {

/// Fill DISTANCE / WAVELENGTH / BEAM / PIXEL from the image header where
/// CONTROL did not set them (marks them ``has_*`` so later steps agree).
void apply_image_geometry_defaults(
    control::ControlState& state,
    const image::DiffractionImage& image
);

/// Image for ``frame_hint``: ``--image`` override, else IMAGE, else TEMPLATE.
/// Throws ``std::runtime_error`` when none is available.
[[nodiscard]] std::filesystem::path resolve_findspots_image_path(
    const control::ControlState& state,
    const std::optional<std::filesystem::path>& image_override,
    int frame_hint
);

/// File holding ``frame`` of a MULTIFRAME wedge. CBF/SMV are one file per
/// frame: the :func:`resolve_findspots_image_path` file is followed to its
/// sibling numbered ``frame`` (last digit run of the stem, same padding),
/// else TEMPLATE; HDF5 / NeXus masters are returned as given. Returning one
/// ``--image`` file for every frame indexed 76c1b from a single 0.2° image.
///
/// :param state: CONTROL state (IMAGE, TEMPLATE, DIRECTORY).
/// :param image_override: ``--image`` path, if any.
/// :param frame: Frame number.
/// :return: Path to read for ``frame``.
/// :raises std::runtime_error: If no file for ``frame`` can be found.
[[nodiscard]] std::filesystem::path resolve_findspots_frame_path(
    const control::ControlState& state,
    const std::optional<std::filesystem::path>& image_override,
    int frame
);

/// Register the DPS AUTOINDEX runner with :class:`control::ControlSession`.
void register_control_autoindex();

/// Register DiffractionLab-style FINDSPOTS (peak find → ``.spt``).
void register_control_findspots();

}  // namespace mosflm::indexing
