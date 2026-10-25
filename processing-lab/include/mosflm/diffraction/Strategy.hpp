#pragma once

#include "mosflm/diffraction/ReflectionPredictor.hpp"

#include <optional>
#include <utility>
#include <vector>

namespace mosflm::diffraction {

/// STRATEGY settings (CONTROL ``STRATEGY``).
struct StrategySettings {
    /// Completeness wanted (fraction): the rotation reported is the smallest
    /// reaching it.
    double target_completeness{0.99};
    /// The target applies to anomalous completeness (Friedel pairs).
    bool anomalous{false};
    /// Image width (°), to turn a rotation into a number of images.
    double osc_deg{0.1};
    /// The φ range collected (°, start < end): rotations are chosen inside
    /// it. Unset: anywhere on a full turn (planning a collection).
    std::optional<std::pair<double, double>> sweep_deg;
    /// φ bin and rotation step (°).
    double step_deg{1.0};
    /// MOSFLM ``STRATEGY AUTO ROTATE n``: also report the best start for a
    /// rotation of ``n``° (0: not asked).
    double rotation_deg{0.0};
    /// Resolution shells in the per-segment statistics.
    int n_shells{8};
};

/// Best completeness for one total rotation.
struct StrategyPoint {
    double rotation_deg{0.0};
    double start_deg{0.0};  ///< Start of the best window.
    double completeness{0.0};  ///< Fraction of unique reflections.
    double anomalous_completeness{0.0};  ///< Fraction of acentric Friedel pairs with both mates.
};

/// One resolution shell of a segment.
struct StrategyShell {
    double d_low_A{0.0};
    double d_high_A{0.0};
    int possible{0};
    int observed{0};
    double completeness{0.0};
    double multiplicity{0.0};
    double anomalous_completeness{0.0};
};

/// One rotation range and what it records.
struct StrategySegment {
    double start_deg{0.0};
    double end_deg{0.0};
    double rotation_deg{0.0};
    int images{0};  ///< ``rotation / osc``, rounded up.
    double completeness{0.0};
    double anomalous_completeness{0.0};
    double multiplicity{0.0};  ///< Observations per observed unique reflection.
    std::vector<StrategyShell> shells;
    /// Frames of the collected sweep (set by the caller, which knows them).
    std::optional<std::pair<int, int>> frames;
};

/// What STRATEGY found.
struct StrategyResult {
    int unique_possible{0};  ///< Unique reflections to the resolution limit.
    int acentric_possible{0};
    double dmin_A{0.0};
    double dmax_A{0.0};
    double osc_deg{0.0};  ///< Image width the image counts use.
    /// Best completeness for each rotation, ``step_deg`` apart, up to a full
    /// turn (or the sweep), stopping once complete.
    std::vector<StrategyPoint> curve;
    /// Smallest rotation reaching ``target_completeness`` (none if no
    /// rotation does).
    std::optional<StrategySegment> target;
    /// The best window for ``rotation_deg``, if asked.
    std::optional<StrategySegment> rotation;
    /// The most any rotation records (the last point of the curve).
    double max_completeness{0.0};
    double max_anomalous_completeness{0.0};
};

/// MOSFLM-style strategy: which rotation range records the most unique data.
///
/// Every reflection to ``geometry.settings.dmin_A`` is followed over a full
/// turn: both Ewald crossings are solved analytically (spindle +X, as
/// :class:`ReflectionPredictor`), and kept when they land on the detector
/// (``use_panel_beam``, margin and ``min_zeta`` as for prediction). Indices
/// are reduced to the reciprocal asymmetric unit with their Friedel sign
/// (MOSFLM ``ROTATE`` / ``COMPLETE``). For each total rotation the window
/// start with the highest completeness is found.
///
/// :param geometry: Crystal (space group), orientation, detector and
///     resolution limits, as for prediction (the scan is ignored).
/// :param settings: Target, image width, collected sweep.
/// :returns: The completeness curve, the target segment and its shells.
[[nodiscard]] StrategyResult compute_strategy(
    const PredictionRequest& geometry,
    const StrategySettings& settings
);

}  // namespace mosflm::diffraction
