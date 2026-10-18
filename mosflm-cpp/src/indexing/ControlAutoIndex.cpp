#include "mosflm/indexing/ControlAutoIndex.hpp"

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/indexing/AutoIndexer.hpp"
#include "mosflm/indexing/SpotList.hpp"

#include <filesystem>
#include <stdexcept>

namespace mosflm::indexing {
namespace {

void run_autoindex(control::ControlState& state) {
    std::filesystem::path spt = state.autoindex.spot_file;
    if (spt.empty()) {
        throw std::runtime_error(
            "CONTROL AUTOINDEX: FILE spot-list path is required"
        );
    }
    if (!spt.is_absolute()) {
        spt = state.directory / spt;
    }
    const SpotList spots = SpotList::read_spt(spt);

    AutoIndexGeometry geo;
    geo.distance_mm = state.distance_mm;
    geo.wavelength_A = state.wavelength_A;
    geo.beam_x_mm = state.beam_x_mm;
    geo.beam_y_mm = state.beam_y_mm;
    geo.pixel_mm = state.pixel_x_mm;
    geo.dmax_A = state.dmax_A;

    AutoIndexSettings settings;
    settings.min_isigi = state.autoindex.min_isigi;
    settings.hkldev = state.autoindex.hkldev;
    settings.nvect = state.autoindex.nvect;
    settings.nvref = state.autoindex.nvref;
    settings.prefer_character = state.autoindex.prefer_character;
    settings.prefer_system = state.autoindex.prefer_system;

    const AutoIndexResult result = AutoIndexer{}.index(spots, geo, settings);
    if (!result.ok || result.ierr != 0) {
        if (result.ierr == -2) {
            throw std::runtime_error(
                "CONTROL AUTOINDEX: only " +
                std::to_string(result.n_spots_used) +
                " spots pass I/σ ≥ " +
                std::to_string(settings.min_isigi) +
                " (need ≥20). Lower AUTOINDEX THRESHOLD or re-run FINDSPOTS "
                "with a lower THRESHOLD / higher MAXSPOTS."
            );
        }
        throw std::runtime_error(
            "CONTROL AUTOINDEX failed (ierr=" + std::to_string(result.ierr) +
            ")"
        );
    }

    state.cell = result.cell;
    state.has_cell = true;
    state.amat = result.amat;
    state.has_matrix = true;
    state.autoindex_system = result.selected_system;
    state.autoindex_character = result.selected_character;
    state.autoindex_n_spots = result.n_spots_used;
}

}  // namespace

/// Called once from ipmosflm / tests to enable CONTROL AUTOINDEX + FINDSPOTS.
void register_control_autoindex() {
    control::ControlSession::register_autoindex_runner(run_autoindex);
    register_control_findspots();
}

}  // namespace mosflm::indexing
