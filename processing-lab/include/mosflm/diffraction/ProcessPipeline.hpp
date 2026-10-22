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
#include <iosfwd>
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
    /// Minimum ``FRACTIONCALC`` to write an MTZ row (DIALS export keeps
    /// partiality ≳ 0.5). Set ≤ 0 to disable the gate.
    double min_fractioncalc{0.5};
    /// Scale ``I``/``σ`` by ``1/FRACTIONCALC`` for partials (full-equivalent).
    bool scale_by_partiality{true};
    /// AIR absorption coefficient (mm⁻¹); 0 disables. Path ≈ detector distance
    /// × sec(θ) through air beyond the crystal–detector baseline.
    double air_absorption_mm_inv{0.0};
    const image::DiffractionImage& image;
};

/// In-memory spot table + integrated observations from one PROCESS call.
/// Average 3D profile of strong reflections, logged like XDS INTEGRATE.LP.
///
/// Background-subtracted counts of strong slices (summation I/σ ≥ 3, FLAG 0)
/// summed on a (φ, row, column) grid around the predicted pixel, per
/// detector region (``kGrid`` × ``kGrid``). φ sections are whole images from
/// the calculated φ, −``kNphi``/2 … +``kNphi``/2.
struct ProfileLog {
    static constexpr int kNphi = 9;
    static constexpr int kGrid = 3;
    int half{0};          ///< Peak half-width (pixels); side = 2·half + 1.
    int width{0};         ///< Image size, for region centres.
    int height{0};
    double osc_deg{0.0};  ///< Image width (φ section size).
    std::vector<double> sum;  ///< [region][φ][row][column]
    std::vector<int> n;       ///< Slices added per region.

    /// Add another log (same ``half``) into this one.
    void add(const ProfileLog& other);
    [[nodiscard]] bool empty() const;
};

/// Print ``log`` like XDS INTEGRATE.LP: spot and rocking σ per detector
/// position, then the average profile scaled to 100.
///
/// The profile is one φ section per block, between ``PROFILE_MAP`` and
/// ``PROFILE_MAP_END``. Each section is ``# phi=<image offset>`` followed by
/// ``side`` lines of ``side`` space-separated integers. A log scanner can
/// split those lines and stack the sections into an array.
void write_profile_log(std::ostream& os, const ProfileLog& log);

struct ProcessResult {
    std::vector<PredictedReflection> predictions;
    std::vector<IntegratedReflection> integrated;
    std::vector<MeasObservation> meas_observations;
    std::vector<ProfileObservation> profile_observations;
    /// This frame's standard-profile learners (for pooling across frames).
    ProfileAccumulator profile_learned;
    /// Simple H/K/L/I/SIGI view (always filled from ``mosflm_reflections``).
    std::vector<crystal::SimpleReflection> mtz_reflections;
    /// Full MOSFLM PROCESS rows (FLAG / XDET / FRACTIONCALC / …).
    std::vector<crystal::MosflmReflection> mosflm_reflections;
    std::size_t n_overlapped{0};
    std::size_t n_profile_ok{0};
    std::size_t n_omitted{0};  ///< IGFLAG 0 / −64 / −128 style skips
    refinement::PostRefinementResult postrefinement{};
    /// Geometry this call integrated with (positional refinement input).
    ScanRange scan{};
    DetectorGeometry detector{};
    double image_beam_x_mm{0.0};
    double image_beam_y_mm{0.0};
    /// Per-image φ range for MTZ batch headers (filled by the frame loop).
    std::vector<crystal::BatchScan> batch_scans;
    /// Strong-reflection 3D profile from this image (rocking mode).
    ProfileLog profile_log{};
};

/// Own the known-indexed pipeline without Fortran GENERATE / REEK / PROCESS.
class ProcessPipeline {
public:
    [[nodiscard]] ProcessResult run(const ProcessRequest& request) const;

    /// Write ``result.mosflm_reflections`` (18-col MOSFLM) via :class:`MtzAdapter`.
    ///
    /// ``wavelength_A`` is stored in MTZ dataset/batch headers (POINTLESS needs
    /// batch records for unmerged input).
    static void write_mtz(
        const std::filesystem::path& path,
        const crystal::CrystalMetadata& crystal,
        const ProcessResult& result,
        double wavelength_A = 0.0
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
