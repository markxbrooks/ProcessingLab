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

/// Inputs for a PREDICT-style vertical slice (stills / single image first).
struct PredictRequest {
    crystal::CrystalMetadata crystal;
    OrientationMatrix orientation;
    DetectorGeometry detector;
    double phi_start_deg{0.0};
    double phi_end_deg{0.0};  ///< Equal to start ⇒ still
    double dmin_A{1.5};
};

/// Skeleton for MOSFLM ``PREDICT``. Algorithm body lands in a later pass;
/// dump helpers match the differential format frozen in FIRST_ROUTINE.md.
class ReflectionPredictor {
public:
    /// Placeholder prediction: empty list until the Fortran kernel is ported.
    [[nodiscard]] std::vector<PredictedReflection> predict(
        const PredictRequest& request
    ) const;

    /// Write oracle-compatible dump for differential tests.
    ///
    /// Format (one reflection per line)::
    ///
    ///   h k l x_mm y_mm phi_beg phi_end full(0|1)
    ///
    /// Header line: ``PREDICT_DUMP n=<count>``.
    static void write_dump(
        const std::filesystem::path& path,
        const std::vector<PredictedReflection>& spots
    );

    /// Parse a dump written by :func:`write_dump` (or an ipmosflm oracle).
    [[nodiscard]] static std::vector<PredictedReflection> read_dump(
        const std::filesystem::path& path
    );
};

}  // namespace mosflm::diffraction
