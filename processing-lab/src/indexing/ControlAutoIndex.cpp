#include "mosflm/indexing/ControlAutoIndex.hpp"

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/image/ImageIO.hpp"
#include "mosflm/indexing/AutoIndexer.hpp"
#include "mosflm/indexing/BravaisSelect.hpp"
#include "mosflm/indexing/SpotList.hpp"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace mosflm::indexing {
namespace {

/// AUTOINDEX from a FILE without FINDSPOTS never saw the image header, so
/// unset geometry would fall back to built-in defaults (DISTANCE 160 mm).
void fill_geometry_from_header(
    control::ControlState& state,
    const std::optional<std::filesystem::path>& image_override
) {
    if (state.has_distance && state.has_wavelength && state.has_beam &&
        state.has_pixel) {
        return;
    }
    const int frame = state.findspots.image_number > 0
                          ? state.findspots.image_number
                          : 1;
    try {
        const auto img_path =
            resolve_findspots_image_path(state, image_override, frame);
        apply_image_geometry_defaults(
            state, image::read_diffraction_image(img_path, frame)
        );
        // Say where geometry came from: the CLI may substitute a demo image.
        std::cout << "  AUTOINDEX: header geometry from " << img_path
                  << ": DISTANCE " << state.distance_mm << " WAVELENGTH "
                  << state.wavelength_A << " BEAM " << state.beam_x_mm << " "
                  << state.beam_y_mm << '\n';
    } catch (const std::exception& ex) {
        std::cout << "  AUTOINDEX: no image header (" << ex.what()
                  << "); using CONTROL/default geometry\n";
    }
    if (!state.has_distance) {
        std::cout << "  AUTOINDEX: WARNING DISTANCE not set; using default "
                  << state.distance_mm << " mm\n";
    }
}

void run_autoindex(
    control::ControlState& state,
    const std::optional<std::filesystem::path>& image_override
) {
    fill_geometry_from_header(state, image_override);
    std::filesystem::path spt = state.autoindex.spot_file;
    if (spt.empty()) {
        throw std::runtime_error(
            "CONTROL AUTOINDEX: FILE spot-list path is required"
        );
    }
    // FINDSPOTS writes to the output directory (``DIRECTORY`` unless one is
    // set); a spot list supplied with the data stays in ``DIRECTORY``.
    if (!spt.is_absolute()) {
        const auto written = control::output_path(state, spt);
        spt = std::filesystem::exists(written) ? written : state.directory / spt;
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
    settings.sepmin = state.autoindex.sepmin;
    settings.cluster_radius = state.autoindex.cluster_radius;
    settings.use_diff_vector = state.autoindex.use_diff_vector;
    settings.use_fft3d = state.autoindex.use_fft3d;
    settings.use_dps = state.autoindex.use_dps;
    settings.use_idxref_refine = state.autoindex.use_idxref_refine;
    settings.min_indexed_frac = state.autoindex.min_indexed_frac;

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
    // PROCESS must predict in the frame these spots were indexed in.
    state.spt_invertx = spots.header.invflag != 0;
    // IDXREF refined beam / distance jointly with AMAT; keep them together.
    if (result.idxref_applied) {
        state.beam_x_mm = result.idxref_beam_x_mm;
        state.beam_y_mm = result.idxref_beam_y_mm;
        state.has_beam = true;
        state.distance_mm = result.idxref_distance_mm;
        state.has_distance = true;
        std::cout << std::fixed << std::setprecision(3)
                  << "  IDXREF: BEAM " << state.beam_x_mm << " "
                  << state.beam_y_mm << "  DISTANCE " << state.distance_mm
                  << " (applied to CONTROL)\n";
    }
    state.autoindex_indexed_frac = result.fit.indexed_frac;
    state.autoindex_rmsd_hkl = result.fit.rmsd_hkl;
    // SUGGLATT → default SYMMETRY when the user did not set SYMMETRY.
    if (!state.has_symmetry && result.selected_symmetry_number > 0) {
        state.symmetry_number = result.selected_symmetry_number;
        state.has_symmetry = true;
    }

    // Mosflm DPSCHOOSE-style lattice-character / penalty table.
    write_bravais_penalty_table(
        std::cout,
        result.bravais,
        result.selected_character,
        result.selected_system,
        result.cell,
        settings.max_penalty
    );
}

}  // namespace

/// Called once from ipmosflm / tests to enable CONTROL AUTOINDEX + FINDSPOTS.
void register_control_autoindex() {
    control::ControlSession::register_autoindex_runner(run_autoindex);
    register_control_findspots();
}

}  // namespace mosflm::indexing
