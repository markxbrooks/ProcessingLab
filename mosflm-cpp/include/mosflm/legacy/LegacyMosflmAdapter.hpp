#pragma once

#include "mosflm/diffraction/ProfileFitter.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/diffraction/SummationIntegrator.hpp"
#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/geometry/UnitCell.hpp"
#include "mosflm/image/DiffractionImage.hpp"

#include <string>
#include <vector>

namespace mosflm::legacy {

/// Narrow façade for future Fortran / C++ interoperability.
///
/// Keep COMMON blocks, column-major layout, and mangling **inside** this
/// layer. Domain code should never see raw ``CELL(6)`` buffers.
class LegacyMosflmAdapter {
public:
    static void to_cell6(const UnitCell& cell, float out[6]);
    static UnitCell from_cell6(const float in[6]);
    [[nodiscard]] static std::string oracle_hint();

    [[nodiscard]] static std::vector<diffraction::PredictedReflection>
    predict_still(const diffraction::PredictRequest& request);

    [[nodiscard]] static std::vector<diffraction::IntegratedReflection>
    integrate_summation(
        const image::DiffractionImage& image,
        const DetectorGeometry& detector,
        const std::vector<diffraction::PredictedReflection>& predictions,
        const diffraction::IntegrateRequest& request = {}
    );

    /// Route 2D Gaussian profile fitting (INTEG3-style stand-in).
    [[nodiscard]] static std::vector<diffraction::ProfileFittedReflection>
    integrate_profile(
        const image::DiffractionImage& image,
        const DetectorGeometry& detector,
        const std::vector<diffraction::PredictedReflection>& predictions,
        const diffraction::ProfileFitRequest& request = {}
    );
};

}  // namespace mosflm::legacy
