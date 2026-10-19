#pragma once

#include "mosflm/diffraction/ProfileFitter.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/diffraction/SummationIntegrator.hpp"
#include "mosflm/geometry/UnitCell.hpp"

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

    [[nodiscard]] static diffraction::PredictionResult predict_still(
        const diffraction::PredictionRequest& request
    );

    [[nodiscard]] static std::vector<diffraction::IntegratedReflection>
    integrate_summation(const diffraction::IntegrationRequest& request);

    /// Route 2D Gaussian profile fitting (INTEG3-style stand-in).
    [[nodiscard]] static std::vector<diffraction::ProfileFittedReflection>
    integrate_profile(const diffraction::ProfileFitRequest& request);
};

}  // namespace mosflm::legacy
