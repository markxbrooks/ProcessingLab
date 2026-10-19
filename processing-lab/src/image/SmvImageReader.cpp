#include "mosflm/image/DiffractionImage.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace mosflm::image {
namespace {

std::string trim(std::string s) {
    while (!s.empty() &&
           (s.front() == ' ' || s.front() == '\t' || s.front() == '\r' ||
            s.front() == '\n')) {
        s.erase(s.begin());
    }
    while (!s.empty() &&
           (s.back() == ' ' || s.back() == '\t' || s.back() == '\r' ||
            s.back() == '\n')) {
        s.pop_back();
    }
    return s;
}

std::unordered_map<std::string, std::string> parse_smv_header(
    const std::string& text
) {
    std::unordered_map<std::string, std::string> kv;
    const auto open = text.find('{');
    const auto close = text.find('}');
    if (open == std::string::npos || close == std::string::npos || close < open) {
        throw std::runtime_error("SmvImageReader: missing {HEADER} braces");
    }
    std::string body = text.substr(open + 1, close - open - 1);
    std::size_t pos = 0;
    while (pos < body.size()) {
        const auto semi = body.find(';', pos);
        if (semi == std::string::npos) {
            break;
        }
        std::string field = trim(body.substr(pos, semi - pos));
        pos = semi + 1;
        if (field.empty()) {
            continue;
        }
        const auto eq = field.find('=');
        if (eq == std::string::npos) {
            continue;
        }
        std::string key = trim(field.substr(0, eq));
        std::string val = trim(field.substr(eq + 1));
        // Upper-case keys for case-insensitive lookup.
        for (char& c : key) {
            c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        }
        kv[key] = val;
    }
    return kv;
}

double require_double(
    const std::unordered_map<std::string, std::string>& kv,
    const char* key
) {
    const auto it = kv.find(key);
    if (it == kv.end()) {
        throw std::runtime_error(
            std::string("SmvImageReader: missing header key ") + key
        );
    }
    try {
        return std::stod(it->second);
    } catch (const std::exception&) {
        throw std::runtime_error(
            std::string("SmvImageReader: bad numeric value for ") + key
        );
    }
}

int require_int(
    const std::unordered_map<std::string, std::string>& kv,
    const char* key
) {
    return static_cast<int>(require_double(kv, key));
}

std::string get_or(
    const std::unordered_map<std::string, std::string>& kv,
    const char* key,
    const char* fallback
) {
    const auto it = kv.find(key);
    return it == kv.end() ? std::string(fallback) : it->second;
}

}  // namespace

DiffractionImage SmvImageReader::read(const std::filesystem::path& path) const {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("SmvImageReader: cannot open " + path.string());
    }

    // Read enough ASCII to locate HEADER_BYTES, then re-seek.
    std::string prefix(4096, '\0');
    in.read(prefix.data(), static_cast<std::streamsize>(prefix.size()));
    prefix.resize(static_cast<std::size_t>(in.gcount()));
    const auto kv = parse_smv_header(prefix);

    const int header_bytes = require_int(kv, "HEADER_BYTES");
    const int size1 = require_int(kv, "SIZE1");
    const int size2 = require_int(kv, "SIZE2");
    if (header_bytes <= 0 || size1 <= 0 || size2 <= 0) {
        throw std::runtime_error("SmvImageReader: invalid SIZE/HEADER_BYTES");
    }

    const std::string type = get_or(kv, "TYPE", "unsigned_short");
    std::string type_lc = type;
    for (char& c : type_lc) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    if (type_lc != "unsigned_short") {
        throw std::runtime_error(
            "SmvImageReader: unsupported TYPE=" + type + " (need unsigned_short)"
        );
    }

    const std::string order = get_or(kv, "BYTE_ORDER", "little_endian");
    std::string order_lc = order;
    for (char& c : order_lc) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    const bool little = order_lc.find("little") != std::string::npos;

    DiffractionImage image;
    image.metadata.width = size1;
    image.metadata.height = size2;
    image.metadata.pixel_size_x_mm = require_double(kv, "PIXEL_SIZE");
    image.metadata.pixel_size_y_mm = image.metadata.pixel_size_x_mm;
    image.metadata.detector_distance_mm = require_double(kv, "DISTANCE");
    image.metadata.wavelength_angstrom = require_double(kv, "WAVELENGTH");
    image.metadata.beam_x_mm = require_double(kv, "BEAM_CENTER_X");
    image.metadata.beam_y_mm = require_double(kv, "BEAM_CENTER_Y");
    image.metadata.oscillation_start_deg =
        kv.count("OSC_START") ? require_double(kv, "OSC_START")
                              : (kv.count("PHI") ? require_double(kv, "PHI") : 0.0);
    image.metadata.oscillation_range_deg =
        kv.count("OSC_RANGE") ? require_double(kv, "OSC_RANGE") : 0.0;
    image.metadata.detector_id =
        kv.count("DETECTOR_SN") ? ("ADSC-" + kv.at("DETECTOR_SN")) : "ADSC-SMV";

    const std::size_t npix =
        static_cast<std::size_t>(size1) * static_cast<std::size_t>(size2);
    std::vector<unsigned char> raw(npix * 2U);
    in.clear();
    in.seekg(header_bytes, std::ios::beg);
    in.read(
        reinterpret_cast<char*>(raw.data()),
        static_cast<std::streamsize>(raw.size())
    );
    if (static_cast<std::size_t>(in.gcount()) != raw.size()) {
        throw std::runtime_error(
            "SmvImageReader: truncated pixel data in " + path.string()
        );
    }

    image.pixels.resize(npix);
    for (std::size_t i = 0; i < npix; ++i) {
        const std::size_t o = i * 2U;
        std::uint16_t v = 0;
        if (little) {
            v = static_cast<std::uint16_t>(raw[o]) |
                (static_cast<std::uint16_t>(raw[o + 1U]) << 8);
        } else {
            v = static_cast<std::uint16_t>(raw[o + 1U]) |
                (static_cast<std::uint16_t>(raw[o]) << 8);
        }
        image.pixels[i] = static_cast<std::int32_t>(v);
    }
    return image;
}

}  // namespace mosflm::image
