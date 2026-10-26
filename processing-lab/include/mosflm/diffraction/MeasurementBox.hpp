#pragma once

#include "mosflm/diffraction/SummationIntegrator.hpp"
#include "mosflm/image/DiffractionImage.hpp"

#include <cmath>
#include <limits>
#include <optional>
#include <utility>
#include <cstdint>
#include <memory>
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
    /// RMS residual of the kept background pixels about the plane
    /// (MOSFLM ``RMSBG``).
    double rms{0.0};

    [[nodiscard]] double at(double dc, double dr) const {
        return a * dc + b * dr + c;
    }
};

/**
 * Highest background count kept for an expected count ``lambda``.
 *
 * Returns the smallest ``k`` with Poisson ``P(X > k) < p``, where ``p`` is
 * the one-sided Gaussian tail beyond ``nsigma``. Sparse photon-counting
 * backgrounds (λ ≪ 1) keep their 1- and 2-count pixels: a symmetric
 * ``nsigma·sd`` clip rejects them and drives the background to zero.
 */
[[nodiscard]] double poisson_upper_limit(double lambda, double nsigma);

/**
 * Fit a tilted background plane in the measurement-box rim.
 *
 * Excludes the peak square and neighbour peak footprints. Optional
 * ``reject_nsigma`` iterates a one-sided (high) outlier clip at the larger
 * of ``reject_nsigma·sd`` and the Poisson limit (BGREJECT stand-in).
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

/// Where the rotation axis crosses the detector, for :class:`ApertureMap`.
struct RotationAxisLine {
    double beam_col{0.0};  ///< Beam position (px).
    double beam_row{0.0};
    bool along_fast{false};  ///< Axis along image columns (Pilatus "X"); else along rows.

    /// Distance (px) of pixel (``col``, ``row``) from the axis line through the beam.
    [[nodiscard]] double distance(double col, double row) const {
        return along_fast ? std::abs(row - beam_row) : std::abs(col - beam_col);
    }
};

/// Fraction of a spot's flux a peak box encloses, by distance from the
/// rotation axis.
///
/// Spots have long, faint tails along the rotation direction (76c1b: a
/// σ ≈ 1 px core with tails to ±10 px), longer the further they are from the
/// axis, which a box sized for precision cuts off. Dividing by the enclosed
/// fraction measured on strong spots (an aperture correction) restores that
/// flux without enlarging the box, which would add background noise.
struct ApertureMap {
    RotationAxisLine axis;
    double max_distance{1.0};  ///< Distance covered by the bins (px).
    int half_box{0};  ///< The peak box the fractions are for.
    /// Enclosed fraction per distance bin (equal widths from 0 to
    /// ``max_distance``); bins with too few strong spots take ``overall``.
    std::vector<double> fraction;
    std::vector<std::size_t> n_spots;
    double overall{1.0};

    /// Fraction at pixel (``col``, ``row``), linear between bin centres.
    [[nodiscard]] double at(double col, double row) const;
};

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
    /// DEFPIX mask (row-major, 1 = ok, ``trusted_width`` wide): a rejected
    /// pixel under the peak rejects the reflection, as a masked one does.
    std::shared_ptr<const std::vector<std::uint8_t>> defpix_mask{};
    /// Aperture correction: intensities are divided by the fraction of the
    /// spot's flux the peak box encloses at that position (unset: none).
    std::optional<ApertureMap> aperture{};
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
    /// :attr:`BackgroundPlane::rms` of the plane above.
    double bg_rms{0.0};
};

/// Intensity-weighted COM refine of predicted pixel; capped shift.
void refine_centroid_com(
    const image::DiffractionImage& image,
    int& col,
    int& row,
    int half,
    double max_shift_px
);

/// Spot extent from :func:`estimate_spot_extent`.
struct SpotExtent {
    int half_x{0};       ///< Peak-region half-width along columns (px).
    int half_y{0};       ///< Peak-region half-width along rows (px).
    double sigma_x{0.0};  ///< Second-moment σ of the average spot (px).
    double sigma_y{0.0};
    std::size_t n_spots{0};
    /// Peak-box half-width to integrate with: ⌈2·mean(σ)⌉ + 1 px.
    int half_box{0};
};

/**
 * Peak-region size from the image, as MOSFLM (measurement-box / raster
 * optimisation from the average strong-spot profile) and XDS (spot size
 * from strong reflections) do, instead of from SEPARATION.
 *
 * Strong spots near ``centres_px`` (column, row) are re-centred on their
 * 3×3 maximum, background-subtracted (window median), normalised and
 * averaged over ±``window_half``. The peak region is the 4-connected area
 * around the centre where that average is ≥ ``threshold_frac`` of its peak
 * (20% ≈ 1.8σ for a Gaussian; on CCP4 insulin this gave ipmosflm's peak
 * size and the best agreement with it, where 3% added background noise);
 * its extent gives ``half_x``/``half_y`` and its second moments
 * ``sigma_x``/``sigma_y``; ``half_box`` = ⌈2·mean(σ)⌉ + 1 px is the box to
 * integrate with. Pixels nearer another predicted spot are masked out.
 * Neighbours sit at varying offsets, so they blur into the average's
 * background rather than joining the region.
 *
 * :param image: Frame to measure.
 * :param centres_px: Predicted spot centres (column, row).
 * :param window_half: Half-size of the averaging window (px).
 * :param min_snr: 3×3 peak signal / noise needed to use a spot.
 * :param threshold_frac: Peak-region threshold as a fraction of the peak.
 * :param min_spots: Fewer strong spots → no estimate.
 * :return: The extent, or ``std::nullopt``.
 */
[[nodiscard]] std::optional<SpotExtent> estimate_spot_extent(
    const image::DiffractionImage& image,
    const std::vector<std::pair<double, double>>& centres_px,
    int window_half = 15,
    double min_snr = 10.0,
    double threshold_frac = 0.2,
    std::size_t min_spots = 20
);

/// :func:`estimate_spot_extent` over several frames, one at a time (one
/// 0.15° Eiger frame has too few strong spots away from the beam, and
/// holding many frames costs ~70 MB each).
class SpotExtentEstimator {
public:
    explicit SpotExtentEstimator(
        int window_half = 15,
        double min_snr = 10.0,
        double threshold_frac = 0.2,
        std::size_t min_spots = 20
    );
    /// Add one frame's strong spots near ``centres_px`` (column, row). Pixels nearer another predicted
    /// spot (in ``centres_px`` or ``others``) are left out of a spot's
    /// average, so neighbours cannot inflate the box.
    void add(
        const image::DiffractionImage& image,
        const std::vector<std::pair<double, double>>& centres_px,
        const std::vector<std::pair<double, double>>* others = nullptr
    );
    /// Extent so far, or ``std::nullopt`` if too few strong spots.
    [[nodiscard]] std::optional<SpotExtent> result() const;
    /// Bin the spots by distance from this axis line (before :meth:`add`), for :meth:`aperture`.
    void set_axis(const RotationAxisLine& axis, int width, int height);
    /// Fraction of the flux a ±``half_box`` box encloses, by distance from
    /// the rotation axis; the total is the average spot summed to
    /// ``window_half − 3`` px (the ring beyond is its background).
    /// ``std::nullopt`` if too few strong spots or no axis was set.
    [[nodiscard]] std::optional<ApertureMap> aperture(int half_box) const;

private:
    struct Acc {
        std::vector<double> avg;  ///< Σ normalised values per cell.
        std::vector<double> cnt;  ///< Spots that owned each cell.
        std::size_t n{0};
    };
    [[nodiscard]] std::optional<SpotExtent> extent_of(const Acc& acc, std::size_t min_spots) const;
    int R_;
    double min_snr_;
    double threshold_frac_;
    std::size_t min_spots_;
    Acc all_;
    std::vector<Acc> bins_;  ///< By distance from the rotation axis.
    std::optional<RotationAxisLine> axis_;
    double max_distance_{1.0};
    [[nodiscard]] double enclosed_fraction(const Acc& acc, int half_box) const;
};

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
