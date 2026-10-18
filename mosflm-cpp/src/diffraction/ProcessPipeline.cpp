#include "mosflm/diffraction/ProcessPipeline.hpp"

namespace mosflm::diffraction {

ProcessResult ProcessPipeline::run(const ProcessRequest& request) const {
    ProcessResult out;

    PredictionRequest pred_req;
    pred_req.crystal = request.crystal;
    pred_req.orientation = request.orientation;
    pred_req.detector = request.detector;
    pred_req.scan = request.scan;
    pred_req.settings = request.prediction;
    out.predictions = ReflectionPredictor{}.predict(pred_req).reflections;

    IntegrationRequest ireq{
        request.image,
        request.detector,
        out.predictions,
        request.integration
    };
    out.integrated = SummationIntegrator{}.integrate(ireq);

    out.mtz_reflections.reserve(out.integrated.size());
    for (const auto& r : out.integrated) {
        if (!r.ok) {
            continue;
        }
        out.mtz_reflections.push_back(
            {r.h, r.k, r.l, static_cast<float>(r.intensity),
             static_cast<float>(r.sigma)}
        );
    }
    return out;
}

void ProcessPipeline::write_mtz(
    const std::filesystem::path& path,
    const crystal::CrystalMetadata& crystal,
    const ProcessResult& result
) {
    crystal::MtzAdapter::write_simple(path, crystal, result.mtz_reflections);
}

ProcessRequest ProcessPipeline::known_indexed_insulin(
    const image::DiffractionImage& image
) {
    Eigen::Matrix3d A;
    A << -0.012237, -0.000358, -0.002626, 0.002191, 0.005616, -0.010974,
        0.001492, -0.011185, -0.005426;

    DetectorGeometry det = image.metadata.to_detector_geometry();
    // Predictor uses beam-relative mm.
    det.beam_x_mm = 0.0;
    det.beam_y_mm = 0.0;
    if (det.distance_mm <= 0.0) {
        det.distance_mm = 160.0;
    }
    if (det.wavelength_A <= 0.0) {
        det.wavelength_A = 0.9790;
    }
    if (det.pixel_size_x_mm <= 0.0) {
        det.pixel_size_x_mm = 0.0816;
        det.pixel_size_y_mm = 0.0816;
    }

    ScanRange scan;
    scan.phi_start_deg = image.metadata.oscillation_start_deg;
    scan.phi_end_deg = image.metadata.oscillation_start_deg +
                       image.metadata.oscillation_range_deg;
    if (scan.phi_end_deg <= scan.phi_start_deg) {
        scan.phi_start_deg = 0.0;
        scan.phi_end_deg = 1.0;
    }

    PredictionSettings prediction;
    prediction.dmin_A = 1.45;
    prediction.ewald_tolerance_Ainv = 0.0025;
    prediction.mosaic_spread_deg = 0.30;
    prediction.phi_step_deg = 0.025;
    prediction.apply_systematic_absences = true;

    IntegrationSettings integration;
    integration.box_half_px = 3;
    integration.bg_frame_px = 2;
    integration.apply_lp = true;
    integration.mosflm_detector_axes = true;
    integration.image_beam_x_mm = image.metadata.beam_x_mm;
    integration.image_beam_y_mm = image.metadata.beam_y_mm;

    return ProcessRequest{
        crystal::CrystalMetadata::from_parameters(
            78.190, 78.190, 78.190, 90.0, 90.0, 90.0, "I 2 3"
        ),
        OrientationMatrix::from_mosflm_orientation_matrix(A),
        det,
        scan,
        prediction,
        integration,
        image
    };
}

}  // namespace mosflm::diffraction
