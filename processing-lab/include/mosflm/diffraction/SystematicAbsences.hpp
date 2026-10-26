#pragma once

#include "mosflm/crystal/MtzAdapter.hpp"

#include <array>
#include <vector>

namespace mosflm::diffraction {

struct AbsenceZoneScore {
    char axis{'a'};
    int period{0};
    double peak_height{0.0};
    double sd{0.0};
    double probability{0.0};
    int observations{0};
};

struct AbsenceAnalysis {
    std::array<AbsenceZoneScore, 3> axes{};
    int suggested_space_group{0};
    int applied_space_group{0};
    double suggestion_confidence{0.0};
    double twofold_a_cc{0.0};
    double twofold_b_cc{0.0};
    int twofold_a_pairs{0};
    int twofold_b_pairs{0};
    bool laue_group_checked{false};
};

/// Score tetragonal twofold intensity relations and Fourier-screen axial
/// I/sigma observations for periodic absences. Returns a screw-axis hypothesis.
[[nodiscard]] AbsenceAnalysis analyze_systematic_absences(
    const std::vector<crystal::MosflmReflection>& reflections,
    bool tetragonal_metric
);

}  // namespace mosflm::diffraction
