#include "mosflm/image/DiffractionImage.hpp"

#include <cbf.h>
#include <cbf_simple.h>

#include <cstdio>
#include <cstring>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>

namespace mosflm::image {
namespace {

void parse_minicbf_comment_header(const std::string& text, ImageMetadata& meta) {
    // Pilatus-style lines inside _array_data.header_contents, e.g.
    //   # Wavelength 0.873128 A
    //   # Detector_distance 0.263775 m
    //   # Pixel_size 172e-6 m x 172e-6 m
    //   # Beam_xy (752.83, 824.79) pixels
    //   # Start_angle 0.0000 deg.
    //   # Angle_increment 0.1000 deg.
    auto grab = [&](const char* pattern, double& out) {
        std::regex re(pattern, std::regex::icase);
        std::smatch m;
        if (std::regex_search(text, m, re) && m.size() > 1) {
            out = std::stod(m[1].str());
        }
    };

    grab(R"(#\s*Wavelength\s+([0-9.eE+-]+)\s*A)", meta.wavelength_angstrom);

    double distance_m = 0.0;
    grab(R"(#\s*Detector_distance\s+([0-9.eE+-]+)\s*m)", distance_m);
    if (distance_m > 0.0) {
        meta.detector_distance_mm = distance_m * 1000.0;
    }

    {
        std::regex re(
            R"(#\s*Pixel_size\s+([0-9.eE+-]+)\s*m\s*x\s*([0-9.eE+-]+)\s*m)",
            std::regex::icase
        );
        std::smatch m;
        if (std::regex_search(text, m, re) && m.size() > 2) {
            meta.pixel_size_x_mm = std::stod(m[1].str()) * 1000.0;
            meta.pixel_size_y_mm = std::stod(m[2].str()) * 1000.0;
        }
    }

    {
        std::regex re(
            R"(#\s*Beam_xy\s*\(\s*([0-9.eE+-]+)\s*,\s*([0-9.eE+-]+)\s*\)\s*pixels)",
            std::regex::icase
        );
        std::smatch m;
        if (std::regex_search(text, m, re) && m.size() > 2) {
            const double bx_px = std::stod(m[1].str());
            const double by_px = std::stod(m[2].str());
            if (meta.pixel_size_x_mm > 0.0) {
                meta.beam_x_mm = bx_px * meta.pixel_size_x_mm;
                meta.beam_y_mm = by_px * meta.pixel_size_y_mm;
            } else {
                meta.beam_x_mm = bx_px;
                meta.beam_y_mm = by_px;
            }
        }
    }

    grab(R"(#\s*Start_angle\s+([0-9.eE+-]+)\s*deg)", meta.oscillation_start_deg);
    grab(R"(#\s*Angle_increment\s+([0-9.eE+-]+)\s*deg)", meta.oscillation_range_deg);

    {
        std::regex re(R"(#\s*Detector:\s*([^\r\n]+))", std::regex::icase);
        std::smatch m;
        if (std::regex_search(text, m, re) && m.size() > 1) {
            meta.detector_id = m[1].str();
            // Trim trailing whitespace
            while (!meta.detector_id.empty() &&
                   (meta.detector_id.back() == ' ' ||
                    meta.detector_id.back() == '\r')) {
                meta.detector_id.pop_back();
            }
        }
    }

    }

std::string read_ascii_prefix(const std::filesystem::path& path, std::size_t nbytes) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("CbfImageReader: cannot open " + path.string());
    }
    std::string buf(nbytes, '\0');
    in.read(buf.data(), static_cast<std::streamsize>(nbytes));
    buf.resize(static_cast<std::size_t>(in.gcount()));
    return buf;
}

/// Pull ``_array_data.header_contents`` text via CBFlib (Dectris/Eiger style).
std::string cbf_header_contents(cbf_handle cbf) {
    std::string out;
    if (cbf_find_category(cbf, "array_data") != 0 ||
        cbf_find_column(cbf, "header_contents") != 0 ||
        cbf_rewind_row(cbf) != 0) {
        return out;
    }
    const char* value = nullptr;
    if (cbf_get_value(cbf, &value) == 0 && value != nullptr) {
        out = value;
    }
    return out;
}

}  // namespace

DiffractionImage CbfImageReader::read(const std::filesystem::path& path) const {
    FILE* fp = std::fopen(path.c_str(), "rb");
    if (!fp) {
        throw std::runtime_error("CbfImageReader: cannot open " + path.string());
    }

    cbf_handle cbf = nullptr;
    if (cbf_make_handle(&cbf) != 0) {
        std::fclose(fp);
        throw std::runtime_error("CbfImageReader: cbf_make_handle failed");
    }

    const int read_err = cbf_read_file(cbf, fp, MSG_DIGEST);
    // cbf_read_file takes ownership of FILE* on success; on failure we close.
    if (read_err != 0) {
        std::fclose(fp);
        cbf_free_handle(cbf);
        throw std::runtime_error(
            "CbfImageReader: cbf_read_file failed (" + std::to_string(read_err) + ")"
        );
    }

    DiffractionImage image;
    ImageMetadata& meta = image.metadata;

    // Best-effort CBFlib wavelength / detector (often missing on miniCBF).
    double wavelength = 0.0;
    if (cbf_get_wavelength(cbf, &wavelength) == 0) {
        meta.wavelength_angstrom = wavelength;
    }

    cbf_detector detector = nullptr;
    if (cbf_construct_detector(cbf, &detector, 0) == 0) {
        double distance = 0.0;
        if (cbf_get_detector_distance(detector, &distance) == 0) {
            // CBFlib returns metres for miniCBF geometry when present.
            meta.detector_distance_mm =
                (distance > 10.0) ? distance : distance * 1000.0;
        }
        cbf_free_detector(detector);
        detector = nullptr;
    }

    size_t ndim1 = 0;
    size_t ndim2 = 0;
    const int size_err = cbf_get_image_size(cbf, 0, 0, &ndim1, &ndim2);
    if (size_err != 0) {
        cbf_free_handle(cbf);
        throw std::runtime_error(
            "CbfImageReader: cbf_get_image_size failed (" +
            std::to_string(size_err) + ")"
        );
    }

    // CBFlib: ndim1 = second (slow) dim, ndim2 = fastest dim for this fixture.
    // Smoke: ndim1=1679 ndim2=1475. Expose width=fast, height=slow.
    meta.height = static_cast<int>(ndim1);
    meta.width = static_cast<int>(ndim2);

    const std::size_t nelem = static_cast<std::size_t>(ndim1) * static_cast<std::size_t>(ndim2);
    image.pixels.assign(nelem, 0);
    const int img_err = cbf_get_image(
        cbf,
        0,
        0,
        image.pixels.data(),
        sizeof(std::int32_t),
        1,
        ndim1,
        ndim2
    );

    // Prefer the structured ``_array_data.header_contents`` text; fall back to
    // scanning the raw ASCII prefix for real files that inline it up front.
    std::string header_text = cbf_header_contents(cbf);
    cbf_free_handle(cbf);
    if (img_err != 0) {
        throw std::runtime_error(
            "CbfImageReader: cbf_get_image failed (" + std::to_string(img_err) + ")"
        );
    }

    if (header_text.empty()) {
        try {
            header_text = read_ascii_prefix(path, 8192);
        } catch (const std::exception&) {
            // Pixels already decoded; metadata remains best-effort.
        }
    }
    try {
        parse_minicbf_comment_header(header_text, meta);
    } catch (const std::exception&) {
        // Metadata stays best-effort.
    }

    return image;
}

}  // namespace mosflm::image
