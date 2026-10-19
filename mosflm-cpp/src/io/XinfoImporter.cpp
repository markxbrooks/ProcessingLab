#include "mosflm/io/XinfoImporter.hpp"

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
    // Preserve ``#`` (not used as xinfo comment; avoids clobbering paths).
    const auto bang = line.find('!');
    if (bang != std::string::npos) {
        line.resize(bang);
    }
    return line;
}

std::vector<std::string> tokenize(const std::string& line) {
    std::istringstream iss(line);
    std::vector<std::string> toks;
    std::string t;
    while (iss >> t) {
        toks.push_back(t);
    }
    return toks;
}

}  // namespace

ImportedSweepGeometry XinfoImporter::parse_string(std::string_view text) {
    ImportedSweepGeometry out;
    bool in_wavelength = false;
    bool in_sweep = false;
    bool got_sweep = false;
    bool got_wavelength_block = false;

    std::istringstream stream{std::string(text)};
    std::string raw;
    while (std::getline(stream, raw)) {
        const std::string line = strip_comment(raw);
        const auto toks = tokenize(line);
        if (toks.empty()) {
            continue;
        }
        const std::string key = upper(toks[0]);

        if (key == "BEGIN" && toks.size() >= 2) {
            const std::string block = upper(toks[1]);
            if (block == "WAVELENGTH" && !got_wavelength_block) {
                in_wavelength = true;
            } else if (block == "SWEEP" && !got_sweep) {
                in_sweep = true;
            }
            continue;
        }
        if (key == "END" && toks.size() >= 2) {
            const std::string block = upper(toks[1]);
            if (block == "WAVELENGTH" && in_wavelength) {
                in_wavelength = false;
                got_wavelength_block = true;
            } else if (block == "SWEEP" && in_sweep) {
                in_sweep = false;
                got_sweep = true;
            }
            continue;
        }

        if (in_wavelength && key == "WAVELENGTH" && toks.size() >= 2) {
            out.has_wavelength = true;
            out.wavelength_A = std::stod(toks[1]);
            continue;
        }

        if (!in_sweep || got_sweep) {
            continue;
        }

        if (key == "WAVELENGTH" && toks.size() >= 2) {
            // Sweep may reference a named wavelength; ignore the name token.
            continue;
        }
        if (key == "DIRECTORY" && toks.size() >= 2) {
            out.has_directory = true;
            out.directory = toks[1];
            for (std::size_t i = 2; i < toks.size(); ++i) {
                out.directory += " ";
                out.directory += toks[i];
            }
        } else if (key == "IMAGE" && toks.size() >= 2) {
            out.has_image = true;
            out.image_path = toks[1];
            for (std::size_t i = 2; i < toks.size(); ++i) {
                out.image_path += toks[i];
            }
        } else if (key == "START_END" && toks.size() >= 3) {
            out.has_process_range = true;
            out.process_first = std::stoi(toks[1]);
            out.process_last = std::stoi(toks[2]);
            if (out.process_last < out.process_first) {
                std::swap(out.process_last, out.process_first);
            }
            if (out.process_last > out.process_first && out.has_image) {
                out.has_template = true;
                out.image_template =
                    template_from_image_name(out.image_path.filename().string());
            }
        }
    }

    if (!out.has_directory && !out.has_image && !out.has_wavelength &&
        !out.has_process_range) {
        throw std::runtime_error(
            "XinfoImporter: no WAVELENGTH/SWEEP fields found"
        );
    }
    return out;
}

ImportedSweepGeometry XinfoImporter::parse_file(
    const std::filesystem::path& path
) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error(
            "XinfoImporter: cannot open " + path.string()
        );
    }
    std::ostringstream oss;
    oss << in.rdbuf();
    return parse_string(oss.str());
}

}  // namespace mosflm::io
