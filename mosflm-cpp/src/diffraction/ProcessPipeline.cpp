#include "mosflm/diffraction/ProcessPipeline.hpp"

#include "mosflm/diffraction/LpCorrection.hpp"
#include "mosflm/geometry/LabFrame.hpp"

#include <algorithm>
#include <cmath>

namespace mosflm::diffraction {
namespace {

/// F77 IGFLAG / IBAD bit mask (WRMTZ write policy).
constexpr int kIbadBgratio = 1;
constexpr int kIbadPkratio = 2;
constexpr int kIbadNegProfile = 4;
constexpr int kIbadBgGradient = 8;
constexpr int kIbadOverlap = 16;
constexpr int kIbadEdge = 32;
constexpr int kIbadFewBg = 64;
constexpr int kOmitOffDetector = 128;

float encode_bgpkratios(double bgrat, double pkrat) {
    const int ia =
        std::min(99, static_cast<int>(std::lround(bgrat * 10.0))) +
        static_cast<int>(std::lround(pkrat * 100.0)) * 100;
    return static_cast<float>(ia);
}

}  // namespace

ProcessResult ProcessPipeline::run(const ProcessRequest& request) const {
    ProcessResult out;

    crystal::CrystalMetadata crystal = request.crystal;
    OrientationMatrix orientation = request.orientation;
    PredictionSettings prediction = request.prediction;

    if (request.enable_postrefinement) {
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

    // Index predictions by (h,k,l) for φ / partial metadata.
    struct PredMeta {
        double phi_begin{};
        double phi_end{};
        bool full{true};
        double x_mm{};
        double y_mm{};
    };
    std::vector<PredMeta> pred_meta;
    pred_meta.reserve(out.predictions.size());
    for (const auto& p : out.predictions) {
        pred_meta.push_back(
            {p.phi_begin_deg, p.phi_end_deg, p.full, p.x_mm, p.y_mm}
        );
    }

    const double beam_x = request.integration.mosflm_detector_axes
                              ? request.integration.image_beam_x_mm
                              : request.image.metadata.beam_x_mm;
    const double beam_y = request.integration.mosflm_detector_axes
                              ? request.integration.image_beam_y_mm
                              : request.image.metadata.beam_y_mm;
    const double px = request.detector.pixel_size_x_mm > 0.0
                          ? request.detector.pixel_size_x_mm
                          : 0.0816;
    const double py = request.detector.pixel_size_y_mm > 0.0
                          ? request.detector.pixel_size_y_mm
                          : 0.0816;

    const std::size_t n = out.integrated.size();
    out.mosflm_reflections.reserve(n);
    out.mtz_reflections.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
        const auto& r = out.integrated[i];
        const PredMeta& pm =
            (i < pred_meta.size())
                ? pred_meta[i]
                : PredMeta{0.0, 0.0, true, r.x_mm, r.y_mm};

        // Unmeasured / failed box → omit (IGFLAG 0 / −128 stand-in).
        if (!r.ok) {
            ++out.n_omitted;
            continue;
        }

        int ibad = 0;
        const ProfileObservation* po = nullptr;
        if (request.use_profile && i < out.profile_observations.size()) {
            po = &out.profile_observations[i];
            if (po->overlapped) {
                ibad |= kIbadOverlap;
            }
            if (po->pkratio > 3.5) {
                ibad |= kIbadPkratio;
            }
            if (po->profile_ok && po->profile_intensity < 0.0) {
                ibad |= kIbadNegProfile;
            }
            if (po->n_bg_rejected > 8) {
                ibad |= kIbadBgratio;
            }
            if (po->npix_bg < 4) {
                ibad |= kIbadFewBg;
            }
        } else if (i < out.meas_observations.size()) {
            const auto& m = out.meas_observations[i];
            if (m.overlapped) {
                ibad |= kIbadOverlap;
            }
            if (m.n_bg_rejected > 8) {
                ibad |= kIbadBgratio;
            }
            if (m.npix_bg < 4) {
                ibad |= kIbadFewBg;
            }
        }

        // Too-few-bg (−64) is omitted by WRMTZ.
        if ((ibad & kIbadFewBg) != 0 && (ibad & ~kIbadFewBg) == 0) {
            ++out.n_omitted;
            continue;
        }

        crystal::MosflmReflection row;
        row.h = r.h;
        row.k = r.k;
        row.l = r.l;
        row.batch = 1;
        // Default I/SIGI from summation; prefer profile when fit succeeded
        // (DIALS-like: profile intensities correlate better for shared HKL).
        row.intensity = static_cast<float>(r.intensity);
        row.sigma = static_cast<float>(r.sigma);
        if (po != nullptr) {
            row.ipr = static_cast<float>(
                po->profile_ok ? po->profile_intensity : po->summation_intensity
            );
            row.sigipr = static_cast<float>(
                po->profile_ok ? po->profile_sigma : po->summation_sigma
            );
            if (po->profile_ok) {
                row.intensity = row.ipr;
                row.sigma = row.sigipr;
            }
            const double bgrat =
                po->npix_bg > 0
                    ? static_cast<double>(po->n_bg_rejected) /
                          static_cast<double>(po->npix_bg)
                    : 0.0;
            row.bgpkratios = encode_bgpkratios(bgrat, po->pkratio);
        } else {
            row.ipr = row.intensity;
            row.sigipr = row.sigma;
        }

        const double phi0 = pm.phi_begin;
        const double phi1 = pm.phi_end;
        const double width = std::max(0.0, phi1 - phi0);
        row.rot = static_cast<float>(0.5 * (phi0 + phi1));
        row.width = static_cast<float>(width);

        // FRACTIONCALC: fraction of the *intrinsic* rocking width on this
        // image (mosaic+divergence). Mid-frame Bragg hits are treated as
        // full even if Ewald sampling touched an osc endpoint.
        //
        // Important: when mosaic ≳ osc (default MOSAIC was 0.30° on 0.15°
        // frames), margin = 0.5·rocking eats the whole frame and every
        // whole-osc observation collapses to osc/mosaic = 0.5 — emptying
        // "fully recorded" zone views (e.g. hk0). Cap the mid-frame margin
        // to a fraction of the oscillation and treat whole-frame interior
        // hits as full.
        const double osc0 = request.scan.phi_start_deg;
        const double osc1 = request.scan.phi_end_deg;
        const double osc_w = std::max(1e-6, osc1 - osc0);
        const double mosaic = std::max(0.0, request.prediction.mosaic_spread_deg);
        const double div = std::max(
            request.prediction.divh_deg, request.prediction.divv_deg
        );
        const double phi_step =
            std::max(1e-6, request.prediction.phi_step_deg);
        const double rocking = std::max(mosaic + div, phi_step);
        const double mid = 0.5 * (phi0 + phi1);
        const double margin = std::min(0.5 * rocking, 0.35 * osc_w);
        bool treat_full = pm.full;
        if (!treat_full && mid - osc0 > margin && osc1 - mid > margin) {
            treat_full = true;
        }
        const double lo = std::max(phi0, osc0);
        const double hi = std::min(phi1, osc1);
        const double overlap = std::max(0.0, hi - lo);
        // Reflection spans (almost) the whole frame with Bragg mid interior.
        if (!treat_full && overlap >= 0.98 * osc_w &&
            mid - osc0 > 0.05 * osc_w && osc1 - mid > 0.05 * osc_w) {
            treat_full = true;
        }
        if (treat_full) {
            row.fractioncalc = 1.0f;
        } else {
            row.fractioncalc = static_cast<float>(
                std::min(1.0, overlap / rocking)
            );
        }
        row.mpart = treat_full ? 0 : 1;

        if (request.min_fractioncalc > 0.0 &&
            static_cast<double>(row.fractioncalc) < request.min_fractioncalc) {
            ++out.n_omitted;
            continue;
        }

        // Full-equivalent intensities: scale partials by 1/FRACTIONCALC so
        // mean-I over HKL matches DIALS-style full-obs comparisons.
        if (request.scale_by_partiality) {
            const double frac =
                std::max(0.05, static_cast<double>(row.fractioncalc));
            if (frac < 0.999) {
                const double inv = 1.0 / frac;
                row.intensity = static_cast<float>(
                    static_cast<double>(row.intensity) * inv
                );
                row.sigma = static_cast<float>(
                    static_cast<double>(row.sigma) * inv
                );
                row.ipr = static_cast<float>(
                    static_cast<double>(row.ipr) * inv
                );
                row.sigipr = static_cast<float>(
                    static_cast<double>(row.sigipr) * inv
                );
            }
        }

        // CORRECT-lite AIR absorption: path length through air beyond the
        // crystal–detector baseline (XDS AIR= stand-in).
        if (request.air_absorption_mm_inv > 0.0 &&
            request.detector.distance_mm > 0.0) {
            const double dx = pm.x_mm - request.detector.beam_x_mm;
            const double dy = pm.y_mm - request.detector.beam_y_mm;
            const double rad = std::hypot(dx, dy);
            const double path =
                std::hypot(request.detector.distance_mm, rad);
            const double extra =
                std::max(0.0, path - request.detector.distance_mm);
            const double corr =
                std::exp(request.air_absorption_mm_inv * extra);
            row.intensity = static_cast<float>(
                static_cast<double>(row.intensity) * corr
            );
            row.sigma = static_cast<float>(
                static_cast<double>(row.sigma) * corr
            );
            row.ipr = static_cast<float>(
                static_cast<double>(row.ipr) * corr
            );
            row.sigipr = static_cast<float>(
                static_cast<double>(row.sigipr) * corr
            );
        }

        double xdet = 0.0;
        double ydet = 0.0;
        // Same YSCALE-aware map as MeasurementBox / StandardProfile.
        cpp_mm_to_mosflm_xdet_yscale(
            pm.x_mm,
            pm.y_mm,
            beam_x,
            beam_y,
            px,
            py,
            request.detector.yscale,
            xdet,
            ydet
        );
        row.xdet = static_cast<float>(xdet);
        row.ydet = static_cast<float>(ydet);

        row.lp = static_cast<float>(lp_factor(
            request.detector, pm.x_mm, pm.y_mm,
            request.integration.polarization_fraction
        ));

        // FLAG: 0 good; else IBAD bits (Aimless rejects non-zero by default).
        row.flag = ibad;

        // ASUPUT-style ASU packing for output indices.
        crystal::MtzAdapter::reduce_to_asu(
            crystal.spacegroup, row.h, row.k, row.l, row.misym
        );

        out.mosflm_reflections.push_back(row);
        out.mtz_reflections.push_back(
            {row.h, row.k, row.l, row.intensity, row.sigma}
        );
    }

    // Silence unused constant warning if edge bit unused for now.
    (void)kIbadBgGradient;
    (void)kIbadEdge;
    (void)kOmitOffDetector;

    return out;
}

void ProcessPipeline::write_mtz(
    const std::filesystem::path& path,
    const crystal::CrystalMetadata& crystal,
    const ProcessResult& result,
    double wavelength_A
) {
    if (!result.mosflm_reflections.empty()) {
        crystal::MtzAdapter::write_mosflm(
            path, crystal, result.mosflm_reflections, wavelength_A
        );
    } else {
        crystal::MtzAdapter::write_simple(
            path, crystal, result.mtz_reflections
        );
    }
}

ProcessRequest ProcessPipeline::known_indexed_insulin(
    const image::DiffractionImage& image
) {
    Eigen::Matrix3d A;
    A << -0.012237, -0.000358, -0.002626, 0.002191, 0.005616, -0.010974,
        0.001492, -0.011185, -0.005426;

    DetectorGeometry det = image.metadata.to_detector_geometry();
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
    prediction.dmax_A = 200.0;
    prediction.ewald_tolerance_Ainv = 0.0025;
    prediction.mosaic_spread_deg = 0.30;
    prediction.phi_step_deg = 0.025;
    prediction.apply_systematic_absences = true;
    prediction.exclude_ice_rings = false;
    prediction.exclude_cusp = true;
    prediction.cusp_fraction = 0.035;
    prediction.max_phi_width_deg = 5.0;
    prediction.apply_reek_detector_margin = true;
    prediction.detector_margin_px = 24.0;
    prediction.widen_ewald_with_mosaic = true;
    prediction.mosaic_ewald_scale = 0.55;

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
    profile.n_bins_x = 5;
    profile.n_bins_y = 5;

    refinement::PostRefinementSettings postref;
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
        // Oracle differentials keep all fractions; CONTROL default is 0.5.
        0.0,
        true,   // scale_by_partiality
        0.0,    // air_absorption_mm_inv
        image
    };
}

}  // namespace mosflm::diffraction
