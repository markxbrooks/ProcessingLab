#pragma once

#include "mosflm/image/DiffractionImage.hpp"

namespace mosflm::diffraction {

/// In-memory INIT-lite estimates (no BLANK.cbf / GAIN.cbf files).
struct DetectorInitEstimate {
    /// Mean corner pedestal (ADU) — BLANK stand-in.
    double blank_adu{0.0};
    /// Robust variance / mean ratio in background corners — GAIN stand-in.
    /// Ideal Poisson photon-counting → ~1; CCD/PAD often 1–3.
    double gain{1.0};
    /// Robust background σ (ADU) in corners after blank subtraction.
    double bg_sigma_adu{1.0};
    bool ok{false};
};

/// Estimate BLANK + GAIN from the four image corners (NBX/NBY-style boxes).
[[nodiscard]] DetectorInitEstimate estimate_blank_gain(
    const image::DiffractionImage& image,
    int corner_box_px = 32
);

}  // namespace mosflm::diffraction
