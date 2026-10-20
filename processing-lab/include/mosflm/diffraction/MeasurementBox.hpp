#pragma once

#include "mosflm/diffraction/SummationIntegrator.hpp"
#include "mosflm/image/DiffractionImage.hpp"

#include <limits>
#include <cstdint>
#include <vector>

namespace mosflm::diffraction {

/// Planar background ``bg = a·Δcol + b·Δrow + c`` (MOSFLM EVAL / INTEG2).
struct BackgroundPlane {
    double a{0.0};
    double b{0.0};
    double c{0.0};
    bool ok{false};
    int n_used{0};
    int n_rejected{0};

    [[nodiscard]] double at(double dc, double dr) const {
        return a * dc + b * dr + c;
    }
};

/**
 * Fit a tilted background plane in the measurement-box rim.
 *
 * Excludes the peak square and neighbour peak footprints. Optional
 * ``reject_nsigma`` iterates residual sigma-clip (BGREJECT stand-in).
 */
[[nodiscard]] BackgroundPlane fit_background_plane(
    const image::DiffractionImage& image,
    int c0,
    int r0,
    int peak_half,
    int outer,
    const std::vector<int>& neigh_cols,
    const std::vector<int>& neigh_rows,
    double reject_nsigma = 3.0,
    const std::vector<std::uint8_t>* trusted_mask = nullptr,
    int trusted_width = 0
);

/// MOSFLM MEAS-inspired measurement-box / background / overlap controls.
///
/// Mirrors the SEPARATION keyword and a simplified BGREJECT / OVERLAP contract
/// from amalgamation ``MEAS``. Peak/background are concentric squares;
/// overlapped neighbours are excluded from the background rim. Background is
/// a tilted plane (EVAL), not a flat mean.
struct MeasurementBoxSettings {
    /// MOSFLM ``SEPARATION`` in millimetres (scanner X / Y).
    double separation_x_mm{0.57};
    double separation_y_mm{0.57};
    /// If >= 0, overrides the peak half-width derived from SEPARATION.
    int peak_half_px{-1};
    /// Background rim thickness outside the peak square (pixels).
    int bg_rim_px{2};
    /// Iterative background outlier clip on plane residuals; 0 disables.
    double bg_reject_nsigma{3.0};
    /// Spots closer than SEPARATION are flagged overlapped; their peak
    /// footprints are excluded from each other's background.
    bool mark_overlaps{true};
    /// If true, overlapped spots are still integrated (``ok`` may stay true)
    /// but ``IntegratedReflection::overlapped`` is set.
    bool integrate_overlaps{true};
    /// Use tilted ``BackgroundPlane`` (default); false → flat mean of rim.
    bool use_bg_plane{true};
    /// Local COM refine of predicted (col,row) before integrate.
    bool refine_centroid{true};
    /// Half-window (px) for COM refine around the predicted pixel.
    int centroid_half_px{3};
    /// Cap on COM shift from prediction (px).
    double centroid_max_shift_px{2.0};
    /// Optional INIT/DEFPIX-lite trusted mask (row-major, 1=ok). Empty = all.
    std::vector<std::uint8_t> trusted_mask{};
    int trusted_width{0};
};

/// Extra fields filled by :class:`MeasIntegrator`.
struct MeasObservation : IntegratedReflection {
    /// Background-subtracted peak centroid in image pixels (column, row);
    /// NaN unless the spot has summation I/σ ≥ 3. With
    /// ``mosflm_detector_axes`` these equal MOSFLM ``XDET``/``YDET``.
    double xobs_px{std::numeric_limits<double>::quiet_NaN()};
    double yobs_px{std::numeric_limits<double>::quiet_NaN()};
    bool overlapped{false};
    /// Peak pixels that were masked (−1) or overloaded (−2); masked ⇒ !ok.
    int n_masked{0};
    int n_overload{0};
    int n_bg_rejected{0};
    int box_half_px{0};
    double bg_a{0.0};
    double bg_b{0.0};
    double bg_c{0.0};
};

/// Intensity-weighted COM refine of predicted pixel; capped shift.
void refine_centroid_com(
    const image::DiffractionImage& image,
    int& col,
    int& row,
    int half,
    double max_shift_px
);

/// Random-access measurement boxes (MEAS science contract without the
/// Fortran stripe I/O buffer). Feeds PROCESS-style summation / profiles.
class MeasIntegrator {
public:
    [[nodiscard]] std::vector<MeasObservation> integrate(
        const IntegrationRequest& request,
        const MeasurementBoxSettings& box
    ) const;
};

}  // namespace mosflm::diffraction
