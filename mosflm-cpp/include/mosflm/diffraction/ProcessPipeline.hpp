#pragma once

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/crystal/MtzAdapter.hpp"
#include "mosflm/diffraction/MeasurementBox.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/diffraction/StandardProfile.hpp"
#include "mosflm/diffraction/SummationIntegrator.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/image/DiffractionImage.hpp"
#include "mosflm/refinement/PostRefinement.hpp"

#include <filesystem>
#include <vector>

namespace mosflm::diffraction {

/// Inputs for a pure-C++ known-indexed PROCESS (predict → MEAS → MTZ).
///
/// ``image`` is referenced (not copied) and must outlive :meth:`ProcessPipeline::run`.
struct ProcessRequest {
    crystal::CrystalMetadata crystal;
    OrientationMatrix orientation;
    /// Predictor detector (beam-relative; typically ``beam_* = 0``).
    DetectorGeometry detector;
    ScanRange scan;
    PredictionSettings prediction{};
    IntegrationSettings integration{};
    /// MEAS-style measurement box / overlap / background (D2).
    MeasurementBoxSettings measurement{};
    /// Standard-profile (INTEG3/4) options (D3).
    StandardProfileSettings profile{};
    /// If true (default), use :class:`StandardProfileFitter` (MEAS + Rossmann).
    /// If false and ``use_meas``, use :class:`MeasIntegrator` only.
    bool use_profile{true};
    /// If true and not ``use_profile``, use :class:`MeasIntegrator`.
    bool use_meas{true};
    /// Stage 2 post-refinement (D4). Off by default until oracle gates pass
    /// for the calling workflow; see CTest ``insulin_postref``.
    bool enable_postrefinement{false};
    refinement::PostRefinementSettings postrefinement{};
    const image::DiffractionImage& image;
};

/// In-memory spot table + integrated observations from one PROCESS call.
struct ProcessResult {
    std::vector<PredictedReflection> predictions;
    std::vector<IntegratedReflection> integrated;
    std::vector<MeasObservation> meas_observations;
    std::vector<ProfileObservation> profile_observations;
    /// Simple H/K/L/I/SIGI view (always filled from ``mosflm_reflections``).
    std::vector<crystal::SimpleReflection> mtz_reflections;
    /// Full MOSFLM PROCESS rows (FLAG / XDET / FRACTIONCALC / …).
    std::vector<crystal::MosflmReflection> mosflm_reflections;
    std::size_t n_overlapped{0};
    std::size_t n_profile_ok{0};
    std::size_t n_omitted{0};  ///< IGFLAG 0 / −64 / −128 style skips
    refinement::PostRefinementResult postrefinement{};
};

/// Own the known-indexed pipeline without Fortran GENERATE / REEK / PROCESS.
class ProcessPipeline {
public:
    [[nodiscard]] ProcessResult run(const ProcessRequest& request) const;

    /// Write ``result.mosflm_reflections`` (18-col MOSFLM) via :class:`MtzAdapter`.
    static void write_mtz(
        const std::filesystem::path& path,
        const crystal::CrystalMetadata& crystal,
        const ProcessResult& result
    );

    /// Geometry / orientation from
    /// ``testdata/images/known_indexed/geometry.txt`` (insulin I23).
    /// Detector comes from ``image.metadata`` (SMV header); prediction uses
    /// beam-centred mm (``beam_* = 0``). Integration uses MOSFLM detector axes
    /// with the image beam centre.
    [[nodiscard]] static ProcessRequest known_indexed_insulin(
        const image::DiffractionImage& image
    );
};

}  // namespace mosflm::diffraction
