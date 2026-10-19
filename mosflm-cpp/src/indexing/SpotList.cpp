#include "mosflm/indexing/SpotList.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
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
        // Multi-image separator (-99 -99 …) — skip, do not treat as a spot.
        if (s.x_mm < -90.0 && s.y_mm < -90.0) {
            continue;
        }
        out.spots.push_back(s);
    }

    if (out.spots.empty()) {
        throw std::runtime_error("SpotList: no spots in " + path.string());
    }
    return out;
}

void SpotList::write_spt(const std::filesystem::path& path) const {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("SpotList: cannot write " + path.string());
    }

    out << std::fixed;
    out << std::setw(12) << header.nfast << std::setw(12) << header.nslow
        << std::setw(11) << std::setprecision(8) << header.pixel_mm
        << std::setw(12) << std::setprecision(6) << header.yscale
        << std::setw(12) << std::setprecision(6) << header.omega_deg << '\n';
    out << std::setw(12) << header.invflag << std::setw(12) << header.liswung
        << '\n';
    out << std::setw(11) << std::setprecision(5) << header.beam_x_mm
        << std::setw(11) << std::setprecision(5) << header.beam_y_mm << '\n';

    for (const auto& s : spots) {
        const double osc = std::fabs(s.phi_end_deg - s.phi_start_deg);
        const double mid = 0.5 * (s.phi_start_deg + s.phi_end_deg);
        // When start==end (still / MOSFLM insulin style), write that φ twice.
        const double c3 = (osc > 1e-9) ? osc : s.phi_start_deg;
        const double c4 = (osc > 1e-9) ? mid : s.phi_end_deg;
        out << std::setw(11) << std::setprecision(2) << s.x_mm
            << std::setw(10) << std::setprecision(2) << s.y_mm
            << std::setw(9) << std::setprecision(3) << c3
            << std::setw(9) << std::setprecision(3) << c4
            << std::setw(12) << std::setprecision(1) << s.intensity
            << std::setw(10) << std::setprecision(1) << s.sigma << '\n';
    }

    out << std::setw(11) << std::setprecision(2) << -999.0
        << std::setw(10) << std::setprecision(2) << -999.0
        << std::setw(9) << std::setprecision(3) << -999.0
        << std::setw(9) << std::setprecision(3) << -999.0
        << std::setw(12) << std::setprecision(1) << -999.0
        << std::setw(10) << std::setprecision(1) << -999.0 << '\n';
    out << std::setw(7) << header.nfast << std::setw(6) << header.nslow
        << std::setw(10) << std::setprecision(4) << header.pixel_mm
        << std::setw(5) << static_cast<int>(header.yscale) << '\n';
}

}  // namespace mosflm::indexing
