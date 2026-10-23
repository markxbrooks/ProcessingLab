#include "mosflm/io/XinfoImporter.hpp"

#include "mosflm/crystal/SpaceGroupInfo.hpp"

#include <cctype>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
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

namespace {

/// One ``BEGIN SWEEP <name>`` … ``END SWEEP`` block.
struct SweepBlock {
    std::string name;
    std::string wavelength_ref;  ///< ``WAVELENGTH <name>`` inside the sweep.
    ImportedSweepGeometry geometry;
    std::size_t crystal{0};  ///< Index into ``XinfoContents::crystals``.
};

struct XinfoContents {
    std::vector<std::pair<std::string, double>> wavelengths;  ///< name, Å
    /// Per ``BEGIN CRYSTAL`` block: ``USER_CELL`` / ``USER_SPACEGROUP``.
    std::vector<ImportedSweepGeometry> crystals{ImportedSweepGeometry{}};
    std::vector<SweepBlock> sweeps;
};

XinfoContents read_contents(std::string_view text) {
    XinfoContents out;
    std::string wl_name;
    bool in_wavelength = false;
    SweepBlock* sweep = nullptr;
    std::size_t crystal = 0;

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
            const std::string name = toks.size() >= 3 ? toks[2] : std::string();
            if (block == "WAVELENGTH") {
                in_wavelength = true;
                wl_name = name;
            } else if (block == "SWEEP") {
                out.sweeps.push_back({name, {}, {}, crystal});
                sweep = &out.sweeps.back();
            } else if (block == "CRYSTAL") {
                out.crystals.emplace_back();
                crystal = out.crystals.size() - 1;
            }
            continue;
        }
        if (key == "END" && toks.size() >= 2) {
            const std::string block = upper(toks[1]);
            if (block == "WAVELENGTH") {
                in_wavelength = false;
            } else if (block == "SWEEP") {
                sweep = nullptr;
            }
            continue;
        }
        if (in_wavelength && key == "WAVELENGTH" && toks.size() >= 2) {
            out.wavelengths.emplace_back(wl_name, std::stod(toks[1]));
            continue;
        }
        if (sweep == nullptr) {
            // Crystal-level user constraints (xia2 ``USER_CELL``,
            // ``USER_SPACEGROUP``).
            auto& c = out.crystals[crystal];
            if (key == "USER_CELL" && toks.size() >= 7) {
                c.has_cell = true;
                c.cell_a = std::stod(toks[1]);
                c.cell_b = std::stod(toks[2]);
                c.cell_c = std::stod(toks[3]);
                c.cell_alpha = std::stod(toks[4]);
                c.cell_beta = std::stod(toks[5]);
                c.cell_gamma = std::stod(toks[6]);
            } else if (key == "USER_SPACEGROUP" && toks.size() >= 2) {
                std::string symbol = toks[1];
                for (std::size_t i = 2; i < toks.size(); ++i) {
                    symbol += " " + toks[i];
                }
                const auto* sg = crystal::SpaceGroupInfo::find_by_name(symbol.c_str());
                if (sg == nullptr) {
                    throw std::runtime_error("XinfoImporter: unknown USER_SPACEGROUP '" + symbol + "'");
                }
                c.has_symmetry = true;
                c.symmetry_number = sg->number;
            }
            continue;
        }
        auto& g = sweep->geometry;
        if (key == "WAVELENGTH" && toks.size() >= 2) {
            sweep->wavelength_ref = toks[1];
        } else if (key == "DIRECTORY" && toks.size() >= 2) {
            g.has_directory = true;
            g.directory = toks[1];
            for (std::size_t i = 2; i < toks.size(); ++i) {
                g.directory += " ";
                g.directory += toks[i];
            }
        } else if (key == "IMAGE" && toks.size() >= 2) {
            g.has_image = true;
            g.image_path = toks[1];
            for (std::size_t i = 2; i < toks.size(); ++i) {
                g.image_path += toks[i];
            }
        } else if (key == "BEAM" && toks.size() >= 3) {
            g.has_beam = true;
            g.beam_x_mm = std::stod(toks[1]);
            g.beam_y_mm = std::stod(toks[2]);
        } else if (key == "DISTANCE" && toks.size() >= 2) {
            g.has_distance = true;
            g.distance_mm = std::stod(toks[1]);
        } else if (key == "START_END" && toks.size() >= 3) {
            g.has_process_range = true;
            g.process_first = std::stoi(toks[1]);
            g.process_last = std::stoi(toks[2]);
            if (g.process_last < g.process_first) {
                std::swap(g.process_last, g.process_first);
            }
        }
    }
    for (auto& sw : out.sweeps) {
        auto& g = sw.geometry;
        if (g.has_process_range && g.process_last > g.process_first && g.has_image) {
            g.has_template = true;
            g.image_template = template_from_image_name(g.image_path.filename().string());
        }
    }
    return out;
}

ImportedSweepGeometry select_sweep(
    const XinfoContents& c, std::optional<std::string_view> sweep_name
) {
    ImportedSweepGeometry out;
    const SweepBlock* sweep = nullptr;
    if (sweep_name) {
        for (const auto& sw : c.sweeps) {
            if (sw.name == *sweep_name) {
                sweep = &sw;
                break;
            }
        }
        if (sweep == nullptr) {
            std::string names;
            for (const auto& sw : c.sweeps) {
                names += (names.empty() ? "" : ", ") + sw.name;
            }
            throw std::runtime_error(
                "XinfoImporter: no sweep '" + std::string(*sweep_name) +
                "' (sweeps: " + (names.empty() ? "none" : names) + ")"
            );
        }
    } else if (!c.sweeps.empty()) {
        sweep = &c.sweeps.front();
    }
    if (sweep != nullptr) {
        out = c.crystals[sweep->crystal];
        out.merge_from(sweep->geometry);
    }
    // The sweep's named wavelength, else the first wavelength block.
    for (const auto& [name, lambda] : c.wavelengths) {
        if (sweep == nullptr || sweep->wavelength_ref.empty() || name == sweep->wavelength_ref) {
            out.has_wavelength = true;
            out.wavelength_A = lambda;
            break;
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

std::string read_text(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error(
            "XinfoImporter: cannot open " + path.string()
        );
    }
    std::ostringstream oss;
    oss << in.rdbuf();
    return oss.str();
}

}  // namespace

ImportedSweepGeometry XinfoImporter::parse_string(std::string_view text) {
    return select_sweep(read_contents(text), std::nullopt);
}

ImportedSweepGeometry XinfoImporter::parse_string(
    std::string_view text, std::string_view sweep
) {
    return select_sweep(read_contents(text), sweep);
}

std::vector<std::string> XinfoImporter::sweep_names(std::string_view text) {
    std::vector<std::string> names;
    for (const auto& sw : read_contents(text).sweeps) {
        names.push_back(sw.name);
    }
    return names;
}

ImportedSweepGeometry XinfoImporter::parse_file(
    const std::filesystem::path& path
) {
    return parse_string(read_text(path));
}

ImportedSweepGeometry XinfoImporter::parse_file(
    const std::filesystem::path& path, std::string_view sweep
) {
    return parse_string(read_text(path), sweep);
}

std::vector<std::string> XinfoImporter::sweep_names_in_file(
    const std::filesystem::path& path
) {
    return sweep_names(read_text(path));
}

}  // namespace mosflm::io
