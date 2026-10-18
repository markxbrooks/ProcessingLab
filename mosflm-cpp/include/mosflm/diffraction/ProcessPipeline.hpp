#pragma once

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/crystal/MtzAdapter.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/diffraction/SummationIntegrator.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/image/DiffractionImage.hpp"

#include <filesystem>
#include <vector>

namespace mosflm::diffraction {

/// Inputs for a pure-C++ known-indexed PROCESS (predict → sum → MTZ).
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
    const image::DiffractionImage& image;
};

/// In-memory spot table + integrated observations from one PROCESS call.
struct ProcessResult {
    std::vector<PredictedReflection> predictions;
    std::vector<IntegratedReflection> integrated;
    std::vector<crystal::SimpleReflection> mtz_reflections;
};

/// Own the known-indexed pipeline without Fortran GENERATE / REEK / PROCESS.
class ProcessPipeline {
public:
    [[nodiscard]] ProcessResult run(const ProcessRequest& request) const;

    /// Write ``result.mtz_reflections`` via :class:`MtzAdapter`.
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
