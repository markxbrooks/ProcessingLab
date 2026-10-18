#pragma once

#include "mosflm/diffraction/SummationIntegrator.hpp"

#include <vector>

namespace mosflm::diffraction {

/// MOSFLM MEAS-inspired measurement-box / background / overlap controls.
///
/// Mirrors the SEPARATION keyword and a simplified BGREJECT / OVERLAP contract
/// from amalgamation ``MEAS`` (~119126). Peak/background are concentric
/// squares; overlapped neighbours are excluded from the background rim.
struct MeasurementBoxSettings {
    /// MOSFLM ``SEPARATION`` in millimetres (scanner X / Y).
    double separation_x_mm{0.57};
    double separation_y_mm{0.57};
    /// If >= 0, overrides the peak half-width derived from SEPARATION.
    int peak_half_px{-1};
    /// Background rim thickness outside the peak square (pixels).
    int bg_rim_px{2};
    /// Iterative background outlier clip (≈ BGREJECT); 0 disables.
    double bg_reject_nsigma{3.0};
    /// Spots closer than SEPARATION are flagged overlapped; their peak
    /// footprints are excluded from each other's background.
    bool mark_overlaps{true};
    /// If true, overlapped spots are still integrated (``ok`` may stay true)
    /// but ``IntegratedReflection::overlapped`` is set.
    bool integrate_overlaps{true};
};

/// Extra fields filled by :class:`MeasIntegrator`.
struct MeasObservation : IntegratedReflection {
    bool overlapped{false};
    int n_bg_rejected{0};
    int box_half_px{0};
};

/// Stripe→shoebox measurement (MEAS contract) feeding PROCESS-style summation.
class MeasIntegrator {
public:
    [[nodiscard]] std::vector<MeasObservation> integrate(
        const IntegrationRequest& request,
        const MeasurementBoxSettings& box
    ) const;
};

}  // namespace mosflm::diffraction
