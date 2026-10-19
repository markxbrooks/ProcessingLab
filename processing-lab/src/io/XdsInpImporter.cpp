#include "mosflm/io/XdsInpImporter.hpp"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mosflm::io {

namespace {

std::string upper(std::string s) {
    for (char& c : s) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return s;
}

std::string strip_comment(std::string line) {
    const auto bang = line.find('!');
    if (bang != std::string::npos) {
        line.resize(bang);
    }
    return line;
}

std::string trim(std::string s) {
    while (!s.empty() &&
           std::isspace(static_cast<unsigned char>(s.front()))) {
        s.erase(s.begin());
    }
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) {
        s.pop_back();
    }
    return s;
}

std::vector<std::string> split_ws(const std::string& s) {
    std::istringstream iss(s);
    std::vector<std::string> toks;
    std::string t;
    while (iss >> t) {
        toks.push_back(t);
    }
    return toks;
}

/// Split a line into KEY and the rest after the first ``=``.
bool split_key_value(const std::string& line, std::string& key, std::string& val) {
    const auto eq = line.find('=');
    if (eq == std::string::npos) {
        return false;
    }
    key = trim(line.substr(0, eq));
    val = trim(line.substr(eq + 1));
    return !key.empty();
}

}  // namespace

ImportedSweepGeometry XdsInpImporter::parse_string(std::string_view text) {
    ImportedSweepGeometry out;
    bool has_orgx = false;
    bool has_orgy = false;
    double orgx = 0.0;
    double orgy = 0.0;

    std::istringstream stream{std::string(text)};
    std::string raw;
    while (std::getline(stream, raw)) {
        const std::string line = trim(strip_comment(raw));
        if (line.empty()) {
            continue;
        }
        std::string key_raw;
        std::string val;
        if (!split_key_value(line, key_raw, val)) {
            continue;
        }
        const std::string key = upper(key_raw);
        const auto toks = split_ws(val);

        if (key == "X-RAY_WAVELENGTH" && !toks.empty()) {
            out.has_wavelength = true;
            out.wavelength_A = std::stod(toks[0]);
        } else if (key == "DETECTOR_DISTANCE" && !toks.empty()) {
            out.has_distance = true;
            out.distance_mm = std::stod(toks[0]);
        } else if (key == "QX" && !toks.empty()) {
            out.pixel_x_mm = std::stod(toks[0]);
            if (!out.has_pixel) {
                out.pixel_y_mm = out.pixel_x_mm;
            }
            out.has_pixel = true;
        } else if (key == "QY" && !toks.empty()) {
            out.pixel_y_mm = std::stod(toks[0]);
            if (!out.has_pixel) {
                out.pixel_x_mm = out.pixel_y_mm;
            }
            out.has_pixel = true;
        } else if (key == "ORGX" && !toks.empty()) {
            has_orgx = true;
            orgx = std::stod(toks[0]);
        } else if (key == "ORGY" && !toks.empty()) {
            has_orgy = true;
            orgy = std::stod(toks[0]);
        } else if (key == "UNIT_CELL_CONSTANTS" && toks.size() >= 6) {
            out.has_cell = true;
            out.cell_a = std::stod(toks[0]);
            out.cell_b = std::stod(toks[1]);
            out.cell_c = std::stod(toks[2]);
            out.cell_alpha = std::stod(toks[3]);
            out.cell_beta = std::stod(toks[4]);
            out.cell_gamma = std::stod(toks[5]);
        } else if (key == "SPACE_GROUP_NUMBER" && !toks.empty()) {
            out.has_symmetry = true;
            out.symmetry_number = std::stoi(toks[0]);
        } else if (key == "DATA_RANGE" && toks.size() >= 2) {
            out.has_process_range = true;
            out.process_first = std::stoi(toks[0]);
            out.process_last = std::stoi(toks[1]);
            if (out.process_last < out.process_first) {
                std::swap(out.process_last, out.process_first);
            }
        } else if (key == "OSCILLATION_RANGE" && !toks.empty()) {
            out.has_oscillation = true;
            out.oscillation_deg = std::stod(toks[0]);
        } else if (key == "NAME_TEMPLATE_OF_DATA_FRAMES" && !toks.empty()) {
            const std::filesystem::path tmpl_path(toks[0]);
            out.has_directory = true;
            out.directory = tmpl_path.parent_path();
            out.has_template = true;
            out.image_template =
                xds_question_marks_to_hash(tmpl_path.filename().string());
        }

        // XDS packs several keys on one line, e.g.
        // NX=2463 NY=2527 QX=0.172 QY=0.172
        // ORGX=... ORGY=...
        // Re-scan the raw line for additional KEY=value tokens.
        std::istringstream multi(line);
        std::string token;
        while (multi >> token) {
            std::string mk;
            std::string mv;
            if (!split_key_value(token, mk, mv)) {
                continue;
            }
            const std::string uk = upper(mk);
            if (uk == "QX") {
                out.pixel_x_mm = std::stod(mv);
                if (!out.has_pixel) {
                    out.pixel_y_mm = out.pixel_x_mm;
                }
                out.has_pixel = true;
            } else if (uk == "QY") {
                out.pixel_y_mm = std::stod(mv);
                if (!out.has_pixel) {
                    out.pixel_x_mm = out.pixel_y_mm;
                }
                out.has_pixel = true;
            } else if (uk == "ORGX") {
                has_orgx = true;
                orgx = std::stod(mv);
            } else if (uk == "ORGY") {
                has_orgy = true;
                orgy = std::stod(mv);
            }
        }
    }

    if (has_orgx && has_orgy && out.has_pixel) {
        out.has_beam = true;
        out.beam_x_mm = orgx * out.pixel_x_mm;
        out.beam_y_mm = orgy * out.pixel_y_mm;
    } else if (has_orgx && has_orgy) {
        // Pixel size missing — leave beam unset (would be wrong in mm).
    }

    if (!out.has_wavelength && !out.has_distance && !out.has_cell &&
        !out.has_process_range && !out.has_template) {
        throw std::runtime_error(
            "XdsInpImporter: no recognised XDS.INP fields found"
        );
    }
    return out;
}

ImportedSweepGeometry XdsInpImporter::parse_file(
    const std::filesystem::path& path
) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error(
            "XdsInpImporter: cannot open " + path.string()
        );
    }
    std::ostringstream oss;
    oss << in.rdbuf();
    return parse_string(oss.str());
}

}  // namespace mosflm::io
