#include "mosflm/indexing/SpotList.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace mosflm::indexing {

SpotList SpotList::read_spt(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("SpotList: cannot open " + path.string());
    }

    SpotList out;
    std::string line;

    if (!std::getline(in, line)) {
        throw std::runtime_error("SpotList: empty file");
    }
    {
        std::istringstream iss(line);
        iss >> out.header.nfast >> out.header.nslow >> out.header.pixel_mm >>
            out.header.yscale >> out.header.omega_deg;
        if (!iss) {
            throw std::runtime_error("SpotList: bad detector header line");
        }
    }

    if (!std::getline(in, line)) {
        throw std::runtime_error("SpotList: missing invflag line");
    }
    {
        std::istringstream iss(line);
        iss >> out.header.invflag >> out.header.liswung;
        if (!iss) {
            throw std::runtime_error("SpotList: bad invflag line");
        }
    }

    if (!std::getline(in, line)) {
        throw std::runtime_error("SpotList: missing beam line");
    }
    {
        std::istringstream iss(line);
        iss >> out.header.beam_x_mm >> out.header.beam_y_mm;
        if (!iss) {
            throw std::runtime_error("SpotList: bad beam line");
        }
    }

    while (std::getline(in, line)) {
        std::istringstream iss(line);
        Spot s;
        if (!(iss >> s.x_mm >> s.y_mm >> s.phi_start_deg >> s.phi_end_deg >>
              s.intensity >> s.sigma)) {
            continue;
        }
        if (s.x_mm < -900.0 && s.y_mm < -900.0) {
            break;  // -999 -999 terminator
        }
        out.spots.push_back(s);
    }

    if (out.spots.empty()) {
        throw std::runtime_error("SpotList: no spots in " + path.string());
    }
    return out;
}

}  // namespace mosflm::indexing
