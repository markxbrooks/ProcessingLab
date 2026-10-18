#include "mosflm/diffraction/ProcessPipeline.hpp"

namespace mosflm::diffraction {

ProcessResult ProcessPipeline::run(const ProcessRequest& request) const {
    ProcessResult out;

    crystal::CrystalMetadata crystal = request.crystal;
    OrientationMatrix orientation = request.orientation;
    PredictionSettings prediction = request.prediction;

    if (request.enable_postrefinement) {
        // Build φ observations from a first-pass prediction (mid φ).
        PredictionRequest probe = {
            crystal,
            orientation,
            request.detector,
            request.scan,
            prediction
        };
        const auto probe_spots =
            ReflectionPredictor{}.predict(probe).reflections;
        std::vector<refinement::PostRefObservation> obs;
        obs.reserve(probe_spots.size());
        for (const auto& s : probe_spots) {
            refinement::PostRefObservation o;
            o.h = s.h;
            o.k = s.k;
            o.l = s.l;
            o.phi_obs_deg = 0.5 * (s.phi_begin_deg + s.phi_end_deg);
            o.weight = s.full ? 0.5 : 1.5;
            obs.push_back(o);
        }
        out.postrefinement = refinement::PostRefiner{}.refine(
            crystal,
            orientation,
            request.detector.wavelength_A,
            prediction.mosaic_spread_deg,
            obs,
            request.postrefinement
        );
        if (out.postrefinement.ok) {
            orientation = out.postrefinement.orientation;
            prediction.mosaic_spread_deg =
                out.postrefinement.parameters.mosaic_spread_deg;
            // Keep crystal metadata; isotropic scale is in orientation UB.
        }
    }

    PredictionRequest pred_req;
    pred_req.crystal = crystal;
    pred_req.orientation = orientation;
    pred_req.detector = request.detector;
    pred_req.scan = request.scan;
    pred_req.settings = prediction;
    out.predictions = ReflectionPredictor{}.predict(pred_req).reflections;

    IntegrationRequest ireq{
        request.image,
        request.detector,
        out.predictions,
        request.integration
    };

    if (request.use_profile) {
        StandardProfileSettings ps = request.profile;
        // Keep box settings in sync with ProcessRequest::measurement.
        ps.box = request.measurement;
        out.profile_observations = StandardProfileFitter{}.fit(ireq, ps);
        out.integrated.reserve(out.profile_observations.size());
        out.meas_observations.reserve(out.profile_observations.size());
        for (const auto& p : out.profile_observations) {
            if (p.overlapped) {
                ++out.n_overlapped;
            }
            if (p.profile_ok) {
                ++out.n_profile_ok;
            }
            out.meas_observations.push_back(p);
            out.integrated.push_back(p);
        }
    } else if (request.use_meas) {
        out.meas_observations =
            MeasIntegrator{}.integrate(ireq, request.measurement);
        out.integrated.reserve(out.meas_observations.size());
        for (const auto& m : out.meas_observations) {
            if (m.overlapped) {
                ++out.n_overlapped;
            }
            out.integrated.push_back(m);
        }
    } else {
        out.integrated = SummationIntegrator{}.integrate(ireq);
    }

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
    // Predictor uses beam-centred mm; image-frame beam is in IntegrationSettings.
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
    // Refined fiducials from oracle PROCESS (insulin_1_001).
    det.yscale = 0.9921;
    det.tilt_mrad = -83.0;
    det.twist_mrad = -140.0;
    det.use_rotmat2 = true;

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

    MeasurementBoxSettings measurement;
    measurement.separation_x_mm = 0.57;
    measurement.separation_y_mm = 0.57;
    measurement.bg_rim_px = 2;
    measurement.bg_reject_nsigma = 3.0;
    measurement.mark_overlaps = true;
    measurement.integrate_overlaps = true;

    StandardProfileSettings profile;
    profile.box = measurement;
    profile.learn_min_isigi = 5.0;
    profile.min_learners = 8;
    profile.fallback_sigma_px = 1.5;
    profile.prefer_profile = true;
    // MOSFLM default above 2.5 Å: 25 standard profiles (5×5 tiling).
    profile.n_bins_x = 5;
    profile.n_bins_y = 5;

    refinement::PostRefinementSettings postref;
    // Opt-in only — CTest insulin_postref exercises the kernel directly.
    const bool enable_postref = false;

    return ProcessRequest{
        crystal::CrystalMetadata::from_parameters(
            78.190, 78.190, 78.190, 90.0, 90.0, 90.0, "I 2 3"
        ),
        OrientationMatrix::from_mosflm_orientation_matrix(A),
        det,
        scan,
        prediction,
        integration,
        measurement,
        profile,
        true,
        true,
        enable_postref,
        postref,
        image
    };
}

}  // namespace mosflm::diffraction
