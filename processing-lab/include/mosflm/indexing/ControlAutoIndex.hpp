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

/// Register the DPS AUTOINDEX runner with :class:`control::ControlSession`.
void register_control_autoindex();

/// Register DiffractionLab-style FINDSPOTS (peak find → ``.spt``).
void register_control_findspots();

}  // namespace mosflm::indexing
