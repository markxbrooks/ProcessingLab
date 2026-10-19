#include "mosflm/io/ImportedSweepGeometry.hpp"

#include <cctype>
#include <sstream>

namespace mosflm::io {

namespace {

std::string format_double(double v) {
    std::ostringstream oss;
    oss.precision(10);
    oss << v;
    return oss.str();
}

}  // namespace

std::string xds_question_marks_to_hash(const std::string& name) {
    std::string out = name;
    for (char& c : out) {
        if (c == '?') {
            c = '#';
        }
    }
    return out;
}

std::string template_from_image_name(const std::string& filename) {
    const auto dot = filename.find_last_of('.');
    const std::size_t end =
        (dot == std::string::npos) ? filename.size() : dot;
    if (end == 0) {
        return filename;
    }
    std::size_t i = end;
    while (i > 0 && std::isdigit(static_cast<unsigned char>(filename[i - 1]))) {
        --i;
    }
    if (i == end) {
        return filename;
    }
    const std::size_t n_digits = end - i;
    std::string out = filename.substr(0, i);
    out.append(n_digits, '#');
    if (dot != std::string::npos) {
        out.append(filename.substr(dot));
    }
    return out;
}

void ImportedSweepGeometry::merge_from(const ImportedSweepGeometry& other) {
    if (other.has_wavelength) {
        has_wavelength = true;
        wavelength_A = other.wavelength_A;
    }
    if (other.has_distance) {
        has_distance = true;
        distance_mm = other.distance_mm;
    }
    if (other.has_pixel) {
        has_pixel = true;
        pixel_x_mm = other.pixel_x_mm;
        pixel_y_mm = other.pixel_y_mm;
    }
    if (other.has_beam) {
        has_beam = true;
        beam_x_mm = other.beam_x_mm;
        beam_y_mm = other.beam_y_mm;
    }
    if (other.has_cell) {
        has_cell = true;
        cell_a = other.cell_a;
        cell_b = other.cell_b;
        cell_c = other.cell_c;
        cell_alpha = other.cell_alpha;
        cell_beta = other.cell_beta;
        cell_gamma = other.cell_gamma;
    }
    if (other.has_symmetry) {
        has_symmetry = true;
        symmetry_number = other.symmetry_number;
    }
    if (other.has_process_range) {
        has_process_range = true;
        process_first = other.process_first;
        process_last = other.process_last;
    }
    if (other.has_oscillation) {
        has_oscillation = true;
        oscillation_deg = other.oscillation_deg;
    }
    if (other.has_directory) {
        has_directory = true;
        directory = other.directory;
    }
    if (other.has_image) {
        has_image = true;
        image_path = other.image_path;
    }
    if (other.has_template) {
        has_template = true;
        image_template = other.image_template;
    }
}

std::vector<std::string> ImportedSweepGeometry::to_control_lines() const {
    std::vector<std::string> lines;
    if (has_wavelength) {
        lines.push_back("WAVELENGTH " + format_double(wavelength_A));
    }
    if (has_distance) {
        lines.push_back("DISTANCE " + format_double(distance_mm));
    }
    if (has_pixel) {
        lines.push_back(
            "PIXEL " + format_double(pixel_x_mm) + " " +
            format_double(pixel_y_mm)
        );
    }
    if (has_beam) {
        lines.push_back(
            "BEAM " + format_double(beam_x_mm) + " " + format_double(beam_y_mm)
        );
    }
    if (has_cell) {
        lines.push_back(
            "CELL " + format_double(cell_a) + " " + format_double(cell_b) +
            " " + format_double(cell_c) + " " + format_double(cell_alpha) +
            " " + format_double(cell_beta) + " " + format_double(cell_gamma)
        );
    }
    if (has_symmetry) {
        lines.push_back("SYMMETRY " + std::to_string(symmetry_number));
    }
    if (has_directory) {
        lines.push_back("DIRECTORY " + directory.string());
    }
    // Multi-frame PROCESS needs TEMPLATE; a lone IMAGE path is reused for
    // every frame by ControlSession, so prefer TEMPLATE when range > 1.
    const bool multi =
        has_process_range && process_last > process_first;
    if (multi && has_template) {
        lines.push_back("TEMPLATE " + image_template);
    } else if (multi && has_image) {
        lines.push_back(
            "TEMPLATE " + template_from_image_name(image_path.filename().string())
        );
    } else {
        if (has_template) {
            lines.push_back("TEMPLATE " + image_template);
        }
        if (has_image) {
            lines.push_back("IMAGE " + image_path.string());
        }
    }
    if (has_process_range) {
        std::ostringstream oss;
        oss << "PROCESS " << process_first << " TO " << process_last;
        if (has_oscillation && oscillation_deg > 0.0) {
            oss << " ANGLE " << format_double(oscillation_deg);
        }
        lines.push_back(oss.str());
    }
    return lines;
}

void apply_imported_geometry(
    control::ControlSession& session,
    const ImportedSweepGeometry& geometry
) {
    for (const auto& line : geometry.to_control_lines()) {
        session.parse_line(line);
    }
}

}  // namespace mosflm::io
