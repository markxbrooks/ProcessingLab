#include "mosflm/indexing/ControlAutoIndex.hpp"

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/diffraction/Defpix.hpp"
#include "mosflm/diffraction/DetectorInit.hpp"
#include "mosflm/image/ImageIO.hpp"
#include "mosflm/indexing/MultiFrameSpotFinder.hpp"
#include "mosflm/indexing/SpotFinder.hpp"
#include "mosflm/indexing/SpotList.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

namespace mosflm::indexing {

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

namespace {

/// HDF5 / NeXus masters hold every frame; the reader selects by index.
bool is_multi_frame_container(const std::filesystem::path& p) {
    std::string ext = p.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return ext == ".h5" || ext == ".hdf5" || ext == ".nxs";
}

/// ``p`` with the last digit run of its stem replaced by ``frame`` (same
/// zero padding), e.g. ``x_1_0001.cbf`` → ``x_1_0037.cbf``; empty if the
/// stem has no digits.
std::filesystem::path with_frame_number(const std::filesystem::path& p, int frame) {
    const std::string stem = p.stem().string();
    const auto end = stem.find_last_of("0123456789");
    if (end == std::string::npos) {
        return {};
    }
    auto begin = end;
    while (begin > 0 && std::isdigit(static_cast<unsigned char>(stem[begin - 1]))) {
        --begin;
    }
    const std::size_t width = end - begin + 1;
    std::string num = std::to_string(frame);
    if (num.size() < width) {
        num.insert(0, width - num.size(), '0');
    }
    return p.parent_path() /
           (stem.substr(0, begin) + num + stem.substr(end + 1) + p.extension().string());
}

/// TEMPLATE (``#`` run → zero-padded ``frame``) under DIRECTORY.
std::filesystem::path template_path(const control::ControlState& state, int frame) {
    const std::string& tmpl = state.image_template;
    const auto hash = tmpl.find('#');
    std::string name = tmpl;
    if (hash != std::string::npos) {
        std::size_t n = 0;
        while (hash + n < tmpl.size() && tmpl[hash + n] == '#') {
            ++n;
        }
        std::string num = std::to_string(frame);
        if (num.size() < n) {
            num.insert(0, n - num.size(), '0');
        }
        name = tmpl.substr(0, hash) + num + tmpl.substr(hash + n);
    }
    return state.directory / name;
}

}  // namespace

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
        return template_path(state, frame_hint);
    }
    throw std::runtime_error(
        "CONTROL FINDSPOTS: no IMAGE/TEMPLATE and no --image override"
    );
}

std::filesystem::path resolve_findspots_frame_path(
    const control::ControlState& state,
    const std::optional<std::filesystem::path>& image_override,
    int frame
) {
    const auto file = resolve_findspots_image_path(state, image_override, frame);
    if (is_multi_frame_container(file)) {
        return file;
    }
    const auto sibling = with_frame_number(file, frame);
    std::error_code ec;
    if (!sibling.empty() && std::filesystem::exists(sibling, ec)) {
        return sibling;
    }
    if (!state.image_template.empty()) {
        return template_path(state, frame);
    }
    throw std::runtime_error(
        "CONTROL FINDSPOTS MULTIFRAME: no image for frame " +
        std::to_string(frame) + " (" + sibling.string() + " missing, no TEMPLATE)"
    );
}

namespace {

/// Frames DEFPIX samples before spot finding (shadow and halo only).
constexpr int kFindspotsDefpixFrames = 5;

/// DEFPIX's pixel classes over ``first``…``last`` for the spot search:
/// shaded, untrusted, defective and hot pixels are left out, and the beam
/// halo is judged against its own background, so the search neither takes
/// the beam's flare for spots nor lets it lift the background statistics.
/// Leaves ``settings`` unchanged with fewer than three frames or with DEFPIX
/// off.
void apply_defpix_to_search(
    const control::ControlState& state,
    const std::optional<std::filesystem::path>& image_override,
    const image::DiffractionImage& probe,
    int first,
    int last,
    SpotFindSettings& search
) {
    if (!state.defpix.enabled || last - first + 1 < 3) {
        return;
    }
    auto settings = state.defpix;
    settings.find_defects = false;  // a spot can stay bright over a few frames
    DetectorGeometry det = probe.metadata.to_detector_geometry();
    det.distance_mm = state.distance_mm;
    det.wavelength_A = state.wavelength_A;
    if (state.has_pixel) {
        det.pixel_size_x_mm = state.pixel_x_mm;
        det.pixel_size_y_mm = state.pixel_y_mm;
    }
    double beam_x = state.beam_x_mm;
    double beam_y = state.beam_y_mm;
    if (!state.has_beam) {
        beam_x = 0.5 * probe.metadata.width * det.pixel_size_x_mm;
        beam_y = 0.5 * probe.metadata.height * det.pixel_size_y_mm;
    }
    const auto frames = diffraction::defpix_frames(first, last, std::min(settings.frames, kFindspotsDefpixFrames));
    const auto res = diffraction::run_defpix(
        [&](int f) {
            return f == first ? probe
                              : image::read_diffraction_image(resolve_findspots_frame_path(state, image_override, f), f);
        },
        frames, det, beam_x, beam_y, settings
    );
    auto mask = std::make_shared<const std::vector<std::uint8_t>>(res.spot_search_mask());
    const auto ignored = static_cast<std::size_t>(std::count(mask->begin(), mask->end(), std::uint8_t{0}));
    const std::size_t n_halo = res.count(diffraction::PixelStatus::Halo);
    if (ignored > 0) {
        search.mask = mask;
    }
    if (n_halo > 0) {
        auto local = std::make_shared<std::vector<float>>(res.status.size(), -1.0f);
        for (std::size_t i = 0; i < res.status.size(); ++i) {
            if (res.status[i] == static_cast<std::uint8_t>(diffraction::PixelStatus::Halo)) {
                (*local)[i] = std::max(0.0f, res.background[i]);
            }
        }
        search.local_background = std::move(local);
    }
    std::cout << "  FINDSPOTS: " << ignored << " pixels left out of the search (shaded "
              << res.count(diffraction::PixelStatus::Shaded) << "); " << n_halo << " in the beam halo";
    if (res.halo.found) {
        std::cout << " (within " << std::fixed << std::setprecision(1) << res.halo.radius_px << " px)";
    }
    std::cout << " judged against their own background; DEFPIX on " << frames.size() << " frames\n";
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
    apply_defpix_to_search(state, image_override, probe, first, last, settings);

    SpotList spots;
    if (state.findspots.multiframe && last > first) {
        MultiFrameSpotSettings mfs;
        mfs.per_frame = settings;
        // A reflection barely moves between frames: 2.5 pixels joins its
        // frames without joining neighbours (spots stay ≥ MINDIST apart).
        constexpr double kMergePx = 2.5;
        mfs.merge_xy_mm = state.findspots.merge_xy_mm < 0.0
                              ? kMergePx * (state.has_pixel ? state.pixel_x_mm : probe.metadata.pixel_size_x_mm)
                              : state.findspots.merge_xy_mm;
        mfs.merge_phi_tol_deg = state.findspots.merge_phi_tol_deg;
        mfs.min_frames = state.findspots.min_frames;
        mfs.max_spots = state.findspots.max_spots;
        mfs.nproc = state.nproc;
        // CBF/SMV ignore the frame argument to read_diffraction_image, so each
        // frame needs its own file (HDF5 masters still use f).
        spots = find_spots_multiframe(
            first, last,
            [&](int f) {
                const auto path_f =
                    resolve_findspots_frame_path(state, image_override, f);
                return image::read_diffraction_image(path_f, f);
            },
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
    out = control::output_path(state, out);
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

    // COLSPOT.LP-style summary (XDS COLSPOT diagnostics).
    const int n_images = (last >= first) ? (last - first + 1) : 1;
    std::cout << "  COLSPOT: images=" << n_images << " (" << first << "-"
              << last << ") spots=" << spots.spots.size()
              << " maxspots=" << state.findspots.max_spots
              << " thr=" << state.findspots.threshold_sigma << " → " << out
              << '\n';
}

}  // namespace

void register_control_findspots() {
    control::ControlSession::register_findspots_runner(run_findspots);
}

}  // namespace mosflm::indexing
