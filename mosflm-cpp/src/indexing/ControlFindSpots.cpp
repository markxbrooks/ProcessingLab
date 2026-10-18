#include "mosflm/indexing/ControlAutoIndex.hpp"

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/image/ImageIO.hpp"
#include "mosflm/indexing/SpotFinder.hpp"
#include "mosflm/indexing/SpotList.hpp"

#include <filesystem>
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

void run_findspots(
    control::ControlState& state,
    const std::optional<std::filesystem::path>& image_override
) {
    // Resolve image path (same rules as ControlSession::resolve_image_path).
    std::filesystem::path img_path;
    if (image_override && !image_override->empty()) {
        img_path = *image_override;
    } else if (!state.image_path.empty()) {
        img_path = state.image_path.is_absolute()
                       ? state.image_path
                       : state.directory / state.image_path;
    } else if (!state.image_template.empty()) {
        // Minimal expand: hash template with findspots image number.
        const std::string& tmpl = state.image_template;
        const auto hash = tmpl.find('#');
        std::string name = tmpl;
        if (hash != std::string::npos) {
            std::size_t n = 0;
            while (hash + n < tmpl.size() && tmpl[hash + n] == '#') {
                ++n;
            }
            std::string num = std::to_string(state.findspots.image_number);
            if (num.size() < n) {
                num.insert(0, n - num.size(), '0');
            }
            name = tmpl.substr(0, hash) + num + tmpl.substr(hash + n);
        }
        img_path = state.directory / name;
    } else {
        throw std::runtime_error(
            "CONTROL FINDSPOTS: no IMAGE/TEMPLATE and no --image override"
        );
    }

    const int frame = state.findspots.image_number > 0
                          ? state.findspots.image_number
                          : 1;
    const auto image = image::read_diffraction_image(img_path, frame);
    apply_image_geometry_defaults(state, image);

    SpotFindSettings settings;
    settings.threshold_sigma = state.findspots.threshold_sigma;
    settings.min_distance = state.findspots.min_distance;
    settings.max_spots = state.findspots.max_spots;
    settings.border = state.findspots.border;
    if (state.findspots.has_phi) {
        settings.phi_start_deg = state.findspots.phi_start_deg;
        settings.phi_end_deg = state.findspots.phi_end_deg;
    }

    SpotList spots = SpotFinder{}.find(image, settings);
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
