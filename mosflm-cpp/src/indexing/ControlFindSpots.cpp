#include "mosflm/indexing/ControlAutoIndex.hpp"

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/diffraction/DetectorInit.hpp"
#include "mosflm/image/ImageIO.hpp"
#include "mosflm/indexing/MultiFrameSpotFinder.hpp"
#include "mosflm/indexing/SpotFinder.hpp"
#include "mosflm/indexing/SpotList.hpp"

#include <algorithm>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>

namespace mosflm::indexing {
namespace {

void apply_image_geometry_defaults(
    control::ControlState& state,
    const image::DiffractionImage& image
) {
    const auto& m = image.metadata;
    if (!state.has_distance && m.detector_distance_mm > 0.0) {
        state.distance_mm = m.detector_distance_mm;
        state.has_distance = true;
    }
    if (!state.has_wavelength && m.wavelength_angstrom > 0.0) {
        state.wavelength_A = m.wavelength_angstrom;
        state.has_wavelength = true;
    }
    if (!state.has_beam && (m.beam_x_mm > 0.0 || m.beam_y_mm > 0.0)) {
        state.beam_x_mm = m.beam_x_mm;
        state.beam_y_mm = m.beam_y_mm;
        state.has_beam = true;
    }
    if (!state.has_pixel && m.pixel_size_x_mm > 0.0) {
        state.pixel_x_mm = m.pixel_size_x_mm;
        state.pixel_y_mm =
            (m.pixel_size_y_mm > 0.0) ? m.pixel_size_y_mm : m.pixel_size_x_mm;
        state.has_pixel = true;
    }
}

std::filesystem::path resolve_findspots_image_path(
    const control::ControlState& state,
    const std::optional<std::filesystem::path>& image_override,
    int frame_hint
) {
    if (image_override && !image_override->empty()) {
        return *image_override;
    }
    if (!state.image_path.empty()) {
        return state.image_path.is_absolute()
                   ? state.image_path
                   : state.directory / state.image_path;
    }
    if (!state.image_template.empty()) {
        const std::string& tmpl = state.image_template;
        const auto hash = tmpl.find('#');
        std::string name = tmpl;
        if (hash != std::string::npos) {
            std::size_t n = 0;
            while (hash + n < tmpl.size() && tmpl[hash + n] == '#') {
                ++n;
            }
            std::string num = std::to_string(frame_hint);
            if (num.size() < n) {
                num.insert(0, n - num.size(), '0');
            }
            name = tmpl.substr(0, hash) + num + tmpl.substr(hash + n);
        }
        return state.directory / name;
    }
    throw std::runtime_error(
        "CONTROL FINDSPOTS: no IMAGE/TEMPLATE and no --image override"
    );
}

void run_findspots(
    control::ControlState& state,
    const std::optional<std::filesystem::path>& image_override
) {
    const int frame = state.findspots.image_number > 0
                          ? state.findspots.image_number
                          : 1;
    int first = frame;
    int last = frame;
    if (state.findspots.multiframe) {
        first = state.findspots.range_first > 0 ? state.findspots.range_first
                                               : frame;
        last = state.findspots.range_last > 0 ? state.findspots.range_last
                                             : first;
        if (last < first) {
            std::swap(first, last);
        }
    }

    const auto img_path =
        resolve_findspots_image_path(state, image_override, first);
    const auto probe = image::read_diffraction_image(img_path, first);
    apply_image_geometry_defaults(state, probe);

    SpotFindSettings settings;
    settings.threshold_sigma = state.findspots.threshold_sigma;
    settings.min_distance = state.findspots.min_distance;
    settings.max_spots = state.findspots.max_spots;
    settings.border = state.findspots.border;
    if (state.findspots.has_phi) {
        settings.phi_start_deg = state.findspots.phi_start_deg;
        settings.phi_end_deg = state.findspots.phi_end_deg;
    }
    {
        const auto init =
            diffraction::estimate_blank_gain(probe, /*corner_box_px=*/32);
        if (init.ok) {
            settings.blank_adu = init.blank_adu;
            settings.gain = init.gain;
        }
    }

    SpotList spots;
    if (state.findspots.multiframe && last > first) {
        MultiFrameSpotSettings mfs;
        mfs.per_frame = settings;
        mfs.merge_xy_mm = state.findspots.merge_xy_mm;
        mfs.merge_phi_tol_deg = state.findspots.merge_phi_tol_deg;
        mfs.min_frames = state.findspots.min_frames;
        mfs.max_spots = state.findspots.max_spots;
        spots = find_spots_multiframe(
            first, last,
            [&](int f) { return image::read_diffraction_image(img_path, f); },
            mfs
        );
    } else {
        spots = SpotFinder{}.find(probe, settings);
    }

    // Prefer CONTROL beam/pixel overrides in the .spt header.
    if (state.has_beam) {
        spots.header.beam_x_mm = state.beam_x_mm;
        spots.header.beam_y_mm = state.beam_y_mm;
    }
    if (state.has_pixel) {
        spots.header.pixel_mm = state.pixel_x_mm;
    }

    std::filesystem::path out = state.findspots.spot_file;
    if (out.empty()) {
        out = "spots.spt";
        state.findspots.spot_file = out;
    }
    if (!out.is_absolute()) {
        out = state.directory / out;
    }
    if (out.extension().empty()) {
        out.replace_extension(".spt");
        state.findspots.spot_file = out.filename();
    }

    if (spots.spots.empty()) {
        throw std::runtime_error(
            "CONTROL FINDSPOTS: no spots found on " + img_path.string()
        );
    }
    spots.write_spt(out);
    state.findspots_n_spots = spots.spots.size();
    state.findspots.spot_file = out;
    // Seed AUTOINDEX FILE when not already set.
    if (state.autoindex.spot_file.empty()) {
        state.autoindex.spot_file = out;
    }
}

}  // namespace

void register_control_findspots() {
    control::ControlSession::register_findspots_runner(run_findspots);
}

}  // namespace mosflm::indexing
