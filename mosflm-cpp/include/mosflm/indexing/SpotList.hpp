#pragma once

#include <filesystem>
#include <vector>

namespace mosflm::indexing {

/// One spot from a MOSFLM ``.spt`` list (mm from image origin).
struct Spot {
    double x_mm{};
    double y_mm{};
    double phi_start_deg{};
    double phi_end_deg{};
    double intensity{};
    double sigma{};
};

/// Header / trailer metadata from a MOSFLM spot file.
struct SpotListHeader {
    int nfast{0};
    int nslow{0};
    double pixel_mm{0.0};
    double yscale{1.0};
    double omega_deg{0.0};
    int invflag{0};
    int liswung{0};
    double beam_x_mm{0.0};
    double beam_y_mm{0.0};
};

/// MOSFLM ``.spt`` spot list (FINDSPOTS / WSPOT format).
struct SpotList {
    SpotListHeader header;
    std::vector<Spot> spots;

    /// Read a MOSFLM spot file (header, beam, spots until ``-999``).
    [[nodiscard]] static SpotList read_spt(const std::filesystem::path& path);

    /// Write a MOSFLM spot file (DiffractionLab / WSPOT layout).
    ///
    /// Spot columns 3–4 are oscillation range and mid-φ (degrees), matching
    /// iMosflm ``Spot::write`` / DiffractionLab ``write_merged_spotfile``.
    void write_spt(const std::filesystem::path& path) const;
};

}  // namespace mosflm::indexing
