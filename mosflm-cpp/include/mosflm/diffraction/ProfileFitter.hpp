#pragma once

#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/diffraction/SummationIntegrator.hpp"
#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/image/DiffractionImage.hpp"

#include <span>
#include <vector>

namespace mosflm::diffraction {

/// One profile-fitted observation (2D Gaussian + flat background).
///
/// This is a first C++ stand-in for MOSFLM ``INTEG3``/``INTEG4`` style fitting —
/// not a line port. Intensity is the discrete sum ``A · Σ w`` over the peak box
/// with ``w = exp(−r² / (2 σ²))``.
struct ProfileFittedReflection {
    int h{};
    int k{};
    int l{};
    double x_mm{};
    double y_mm{};
    double intensity{};
    double sigma{};
    double amplitude{};   ///< Peak height A (ADU)
    double background{};  ///< Flat background B (ADU)
    double sigma_px{};    ///< Fitted / fixed Gaussian σ (pixels)
    double chi2{};        ///< Mean squared residual per pixel
    int npix{};
    bool ok{false};
};

/// Algorithm options for 2D Gaussian profile fit.
struct ProfileFitSettings {
    int box_half_px{5};
    /// Fixed isotropic σ in pixels (free-σ fit deferred).
    double sigma_px{1.5};
    /// If true, refine σ by a short 1-D golden search around ``sigma_px``.
    bool refine_sigma{false};
    /// If true, scale intensity/sigma by unpolarized ``lp_factor``.
    bool apply_lp{false};
};

/// Inputs for one profile-fit call (refs must outlive the call).
struct ProfileFitRequest {
    const image::DiffractionImage& image;
    const DetectorGeometry& detector;
    std::span<const PredictedReflection> reflections;
    ProfileFitSettings settings{};
};

/// Box-local 2D Gaussian profile fitter (peak centre fixed at prediction).
class ProfileFitter {
public:
    [[nodiscard]] std::vector<ProfileFittedReflection> fit(
        const ProfileFitRequest& request
    ) const;
};

}  // namespace mosflm::diffraction
