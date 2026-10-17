#pragma once

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace mosflm::diffraction {

/// One predicted spot (I/O contract from docs/FIRST_ROUTINE.md).
struct PredictedReflection {
    int h{};
    int k{};
    int l{};
    double x_mm{};
    double y_mm{};
    double phi_begin_deg{};
    double phi_end_deg{};
    bool full{true};  ///< false ⇒ partial
};

/// Inputs for PREDICT (stills or oscillation).
struct PredictRequest {
    crystal::CrystalMetadata crystal;
    OrientationMatrix orientation;
    DetectorGeometry detector;
    double phi_start_deg{0.0};
    double phi_end_deg{0.0};  ///< Equal to start ⇒ still
    double dmin_A{1.5};
    /// Accept reflections with ``||s1| - 1/λ| ≤`` this (Å⁻¹).
    double ewald_tolerance_Ainv{0.0025};
    /// φ sampling step for oscillation scans (degrees).
    double phi_step_deg{0.05};
};

/// Reflection prediction for stills and fine-sliced oscillations.
///
/// Frame: incident ``s0 = (0,0,1/λ)``, spindle rotation about +X, flat detector
/// at ``z = distance`` via :meth:`DetectorGeometry::lab_to_detector_mm`.
///
/// Oscillation: sample φ in ``[phi_start, phi_end]``; a spot is recorded if it
/// meets the Ewald test on the detector for any sample. ``full`` is true when
/// the on-sphere φ interval lies strictly inside the oscillation (does not
/// touch either endpoint); stills are always ``full``.
class ReflectionPredictor {
public:
    [[nodiscard]] std::vector<PredictedReflection> predict(
        const PredictRequest& request
    ) const;

    static void write_dump(
        const std::filesystem::path& path,
        const std::vector<PredictedReflection>& spots
    );

    [[nodiscard]] static std::vector<PredictedReflection> read_dump(
        const std::filesystem::path& path
    );

    /// Frozen analytical stills request (``predict_still.dump``).
    [[nodiscard]] static PredictRequest analytical_still_fixture();

    /// Frozen analytical oscillation request (``predict_osc.dump``):
    /// same cell/detector as stills, ``φ ∈ [0, 0.5]°``, ``dmin = 5 Å``.
    [[nodiscard]] static PredictRequest analytical_oscillation_fixture();
};

}  // namespace mosflm::diffraction
