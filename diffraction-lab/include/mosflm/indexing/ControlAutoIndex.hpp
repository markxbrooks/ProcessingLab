#pragma once

namespace mosflm::indexing {

/// Register the DPS AUTOINDEX runner with :class:`control::ControlSession`.
void register_control_autoindex();

/// Register DiffractionLab-style FINDSPOTS (peak find → ``.spt``).
void register_control_findspots();

}  // namespace mosflm::indexing
