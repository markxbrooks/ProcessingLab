#include "mosflm/diffraction/ProcessPipeline.hpp"

#include "mosflm/diffraction/LpCorrection.hpp"
#include "mosflm/geometry/LabFrame.hpp"

#include <algorithm>
#include <ostream>
#include <numeric>
#include <iomanip>
#include <array>
#include <cmath>
#include <limits>

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
/// Overloaded (saturated) pixel under the peak.
constexpr int kIbadOverload = 256;

/// Rotation about the predictor spindle (+X), same sense as
/// :class:`ReflectionPredictor`.
Eigen::Matrix3d rotation_about_x_deg(double phi_deg) {
    const double r = phi_deg * (3.14159265358979323846 / 180.0);
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
    R(1, 1) = std::cos(r);
    R(1, 2) = -std::sin(r);
    R(2, 1) = std::sin(r);
    R(2, 2) = std::cos(r);
    return R;
}

float encode_bgpkratios(double bgrat, double pkrat) {
    const int ia =
        std::min(99, static_cast<int>(std::lround(bgrat * 10.0))) +
        static_cast<int>(std::lround(pkrat * 100.0)) * 100;
    return static_cast<float>(ia);
}

/// Add one strong slice to ``log``: background-subtracted pixels around the
/// predicted pixel (``xdet``, ``ydet``) in φ section ``dphi_images``.
void accumulate_profile(
    ProfileLog& log,
    const image::DiffractionImage& image,
    double xdet,
    double ydet,
    int half,
    const MeasObservation& m,
    double dphi_images,
    double osc_deg
) {
    constexpr int kNphi = ProfileLog::kNphi;
    constexpr int kGrid = ProfileLog::kGrid;
    const int ip = static_cast<int>(std::lround(dphi_images)) + kNphi / 2;
    if (half < 1 || ip < 0 || ip >= kNphi) {
        return;
    }
    const int w = image.metadata.width;
    const int h = image.metadata.height;
    const int c0 = static_cast<int>(std::lround(xdet));
    const int r0 = static_cast<int>(std::lround(ydet));
    if (c0 - half < 0 || r0 - half < 0 || c0 + half >= w || r0 + half >= h) {
        return;
    }
    if (log.sum.empty()) {
        const int side = 2 * half + 1;
        log.half = half;
        log.width = w;
        log.height = h;
        log.osc_deg = osc_deg;
        log.sum.assign(
            static_cast<std::size_t>(kGrid * kGrid * kNphi * side * side), 0.0
        );
        log.n.assign(static_cast<std::size_t>(kGrid * kGrid), 0);
    } else if (half != log.half) {
        return;
    }
    const auto px = [&](int c, int r) {
        return image.pixels[static_cast<std::size_t>(r) * static_cast<std::size_t>(w) +
                            static_cast<std::size_t>(c)];
    };
    for (int dr = -half; dr <= half; ++dr) {
        for (int dc = -half; dc <= half; ++dc) {
            if (px(c0 + dc, r0 + dr) < 0) {
                return;  // masked / overloaded pixel in the box
            }
        }
    }
    const int side = 2 * half + 1;
    const int reg = std::min(kGrid - 1, c0 * kGrid / w) +
                    kGrid * std::min(kGrid - 1, r0 * kGrid / h);
    double* dst = &log.sum[static_cast<std::size_t>((reg * kNphi + ip) * side * side)];
    for (int dr = -half; dr <= half; ++dr) {
        for (int dc = -half; dc <= half; ++dc) {
            const double b = m.bg_a * dc + m.bg_b * dr + m.bg_c;
            dst[(dr + half) * side + (dc + half)] +=
                static_cast<double>(px(c0 + dc, r0 + dr)) - b;
        }
    }
    ++log.n[static_cast<std::size_t>(reg)];
}

}  // namespace

void ProfileLog::add(const ProfileLog& other) {
    if (other.sum.empty()) {
        return;
    }
    if (sum.empty()) {
        *this = other;
        return;
    }
    if (other.half != half) {
        return;
    }
    for (std::size_t i = 0; i < sum.size(); ++i) {
        sum[i] += other.sum[i];
    }
    for (std::size_t i = 0; i < n.size(); ++i) {
        n[i] += other.n[i];
    }
}

bool ProfileLog::empty() const {
    return sum.empty() || std::accumulate(n.begin(), n.end(), 0) == 0;
}

void write_profile_log(std::ostream& os, const ProfileLog& log) {
    if (log.empty()) {
        return;
    }
    constexpr int kNphi = ProfileLog::kNphi;
    constexpr int kGrid = ProfileLog::kGrid;
    constexpr int kRegions = kGrid * kGrid;
    const int side = 2 * log.half + 1;
    const int np = side * side;
    const auto at = [&](int reg, int ip, int k) {
        return log.sum[static_cast<std::size_t>((reg * kNphi + ip) * np + k)];
    };
    const int n_total = std::accumulate(log.n.begin(), log.n.end(), 0);

    // Spot σ (pixels) and rocking σ (degrees) per region from the moments of
    // the positive part of the summed profile (φ binned in whole images).
    std::array<double, kRegions> sxy{};
    std::array<double, kRegions> sphi{};
    for (int reg = 0; reg < kRegions; ++reg) {
        double w = 0.0, mx = 0.0, my = 0.0, mp = 0.0;
        double vx = 0.0, vy = 0.0, vp = 0.0;
        for (int pass = 0; pass < 2; ++pass) {
            for (int ip = 0; ip < kNphi; ++ip) {
                for (int k = 0; k < np; ++k) {
                    const double v = std::max(0.0, at(reg, ip, k));
                    const double x = k % side - log.half;
                    const double y = k / side - log.half;
                    const double f = ip - kNphi / 2;
                    if (pass == 0) {
                        w += v;
                        mx += v * x;
                        my += v * y;
                        mp += v * f;
                    } else {
                        vx += v * (x - mx) * (x - mx);
                        vy += v * (y - my) * (y - my);
                        vp += v * (f - mp) * (f - mp);
                    }
                }
            }
            if (pass == 0) {
                if (!(w > 0.0)) {
                    break;
                }
                mx /= w;
                my /= w;
                mp /= w;
            }
        }
        if (w > 0.0) {
            sxy[static_cast<std::size_t>(reg)] = std::sqrt(0.5 * (vx + vy) / w);
            sphi[static_cast<std::size_t>(reg)] = std::sqrt(vp / w) * log.osc_deg;
        }
    }

    os << "\n SPOT SIZE AND ROCKING WIDTH FROM " << std::setw(7) << n_total
       << " STRONG REFLECTION SLICES AT " << kRegions
       << " POSITIONS ON THE DETECTOR SURFACE.\n"
       << " (SIGMA_XY in pixels; SIGMA_PHI in degrees, phi binned in whole"
          " images of " << std::fixed << std::setprecision(3) << log.osc_deg
       << " deg)\n";
    os << " POSITION NUMBER     ";
    for (int reg = 0; reg < kRegions; ++reg) {
        os << std::setw(7) << reg + 1;
    }
    os << "\n X-COORDINATE (pixel)" << std::setprecision(1);
    for (int reg = 0; reg < kRegions; ++reg) {
        os << std::setw(7) << (reg % kGrid + 0.5) * log.width / kGrid;
    }
    os << "\n Y-COORDINATE (pixel)";
    for (int reg = 0; reg < kRegions; ++reg) {
        os << std::setw(7) << (reg / kGrid + 0.5) * log.height / kGrid;
    }
    os << "\n NUMBER              ";
    for (int reg = 0; reg < kRegions; ++reg) {
        os << std::setw(7) << log.n[static_cast<std::size_t>(reg)];
    }
    os << "\n SIGMA_XY (pixel)    " << std::setprecision(2);
    for (int reg = 0; reg < kRegions; ++reg) {
        os << std::setw(7) << sxy[static_cast<std::size_t>(reg)];
    }
    os << "\n SIGMA_PHI (degree)  " << std::setprecision(3);
    for (int reg = 0; reg < kRegions; ++reg) {
        os << std::setw(7) << sphi[static_cast<std::size_t>(reg)];
    }
    os << "\n";

    // Average profile over all regions, scaled to a maximum of 100.
    std::vector<double> avg(static_cast<std::size_t>(kNphi * np), 0.0);
    for (int reg = 0; reg < kRegions; ++reg) {
        for (int i = 0; i < kNphi * np; ++i) {
            avg[static_cast<std::size_t>(i)] += at(reg, i / np, i % np);
        }
    }
    const double peak = *std::max_element(avg.begin(), avg.end());
    // One φ section per block. Width-4 fields stay separated at the scale
    // peak (100); the old width-3 dump glued 41 and 100 into "41100".
    // PROFILE_MAP / PROFILE_MAP_END bracket the block for a log scanner.
    os << "\n\n ***** AVERAGE THREE-DIMENSIONAL PROFILE OF " << std::setw(7)
       << n_total << " STRONG REFLECTION SLICES *****\n"
       << " (one section per image, " << -kNphi / 2 << " ... +" << kNphi / 2
       << " from the calculated phi; " << side << "x" << side
       << " pixels around the predicted position; scaled to 100)\n"
       << "PROFILE_MAP side=" << side << " nphi=" << kNphi << " scale=100\n";
    for (int ip = 0; ip < kNphi; ++ip) {
        os << "# phi=" << (ip - kNphi / 2) << '\n';
        for (int r = 0; r < side; ++r) {
            for (int c = 0; c < side; ++c) {
                const double v =
                    avg[static_cast<std::size_t>(ip * np + r * side + c)];
                const int scaled =
                    peak > 0.0
                        ? static_cast<int>(std::lround(100.0 * v / peak))
                        : 0;
                if (c != 0) {
                    os << ' ';
                }
                os << std::setw(4) << scaled;
            }
            os << '\n';
        }
    }
    os << "PROFILE_MAP_END\n";
    os << std::flush;
}

ProcessResult ProcessPipeline::run(const ProcessRequest& request) const {
    ProcessResult out;
    out.scan = request.scan;
    out.detector = request.detector;
    out.image_beam_x_mm = request.integration.image_beam_x_mm;
    out.image_beam_y_mm = request.integration.image_beam_y_mm;

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
        out.profile_observations =
            StandardProfileFitter{}.fit(ireq, ps, &out.profile_learned);
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
        double phi_calc{std::numeric_limits<double>::quiet_NaN()};
        double fraction{1.0};
    };
    std::vector<PredMeta> pred_meta;
    pred_meta.reserve(out.predictions.size());
    for (const auto& p : out.predictions) {
        pred_meta.push_back(
            {p.phi_begin_deg, p.phi_end_deg, p.full, p.x_mm, p.y_mm,
             p.phi_calc_deg, p.fraction}
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
            // Non-flat background: many rim pixels rejected. Relative to the
            // rim size (Eiger rims hold ~240 px, where 8 is only ~3%).
            if (po->n_bg_rejected > std::max(8, po->npix_bg / 10)) {
                ibad |= kIbadBgratio;
            }
            if (po->npix_bg < 4) {
                ibad |= kIbadFewBg;
            }
            if (po->n_overload > 0) {
                ibad |= kIbadOverload;
            }
        } else if (i < out.meas_observations.size()) {
            const auto& m = out.meas_observations[i];
            if (m.overlapped) {
                ibad |= kIbadOverlap;
            }
            if (m.n_bg_rejected > std::max(8, m.npix_bg / 10)) {
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
        row.h0 = r.h;
        row.k0 = r.k;
        row.l0 = r.l;
        row.batch = 1;
        // MOSFLM convention: I/SIGI = summation, IPR/SIGIPR = profile fit
        // (summation when the fit failed). The per-slice 2D profile fit
        // scatters ~15% vs pixel sums on fine-sliced Eiger data, while
        // summation agrees to ~2%, so I must not be overwritten by IPR.
        row.intensity = static_cast<float>(r.intensity);
        row.sigma = static_cast<float>(r.sigma);
        if (po != nullptr) {
            row.intensity = static_cast<float>(po->summation_intensity);
            row.sigma = static_cast<float>(po->summation_sigma);
            row.ipr = static_cast<float>(
                po->profile_ok ? po->profile_intensity : po->summation_intensity
            );
            row.sigipr = static_cast<float>(
                po->profile_ok ? po->profile_sigma : po->summation_sigma
            );
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

        const bool rocking = request.prediction.rocking_model &&
                             std::isfinite(pm.phi_calc);
        // φ at which the reflection is in diffracting position.
        const double mid =
            rocking ? pm.phi_calc : 0.5 * (pm.phi_begin + pm.phi_end);
        if (rocking) {
            // One slice of a rocking curve: summed across frames later
            // (sum_partial_slices), so no 1/FRACTIONCALC scaling here.
            row.rot = static_cast<float>(pm.phi_calc);
            row.width = static_cast<float>(pm.phi_end - pm.phi_begin);
            row.fractioncalc = static_cast<float>(pm.fraction);
            row.mpart = pm.fraction < 0.99 ? 1 : 0;
            row.phi_obs = static_cast<float>(
                0.5 * (request.scan.phi_start_deg + request.scan.phi_end_deg)
            );
        } else {
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
        // XDS-style average 3D profile of strong slices (logged per batch).
        if (request.integration.mosflm_detector_axes && rocking &&
            po != nullptr && ibad == 0 && po->summation_sigma > 0.0 &&
            po->summation_intensity >= 3.0 * po->summation_sigma) {
            const double osc =
                request.scan.phi_end_deg - request.scan.phi_start_deg;
            if (osc > 0.0) {
                accumulate_profile(
                    out.profile_log, request.image, xdet, ydet,
                    po->box_half_px, *po,
                    (0.5 * (request.scan.phi_start_deg + request.scan.phi_end_deg) -
                     pm.phi_calc) / osc,
                    osc
                );
            }
        }

        // Observed centroid: pixel (col,row) == XDET/YDET on MOSFLM axes.
        if (request.integration.mosflm_detector_axes && po != nullptr &&
            std::isfinite(po->xobs_px) && std::isfinite(po->yobs_px)) {
            row.xdet_obs = static_cast<float>(po->xobs_px);
            row.ydet_obs = static_cast<float>(po->yobs_px);
        }

        // Aperture correction: the peak box holds only this fraction of a
        // spot with long tails (measured on strong spots at this distance
        // from the rotation axis).
        if (request.measurement.aperture) {
            const double f = request.measurement.aperture->at(xdet, ydet);
            if (f > 0.0 && f < 1.0) {
                row.intensity = static_cast<float>(row.intensity / f);
                row.sigma = static_cast<float>(row.sigma / f);
                row.ipr = static_cast<float>(row.ipr / f);
                row.sigipr = static_cast<float>(row.sigipr / f);
            }
        }

        if (request.integration.rotation_lp) {
            // Rotation-method LP: s1 from UB·h rotated to the reflecting φ.
            const double lam = request.detector.wavelength_A;
            const Eigen::Vector3d s0(0.0, 0.0, lam > 0.0 ? 1.0 / lam : 1.0);
            const Eigen::Vector3d s1 =
                s0 + rotation_about_x_deg(mid) * orientation.miller_to_lab(r.h, r.k, r.l);
            const double corr = rotation_lp_correction(
                s0, s1, Eigen::Vector3d::UnitX(), Eigen::Vector3d::UnitY(),
                request.integration.polarization_fraction
            );
            row.intensity = static_cast<float>(row.intensity * corr);
            row.sigma = static_cast<float>(row.sigma * corr);
            row.ipr = static_cast<float>(row.ipr * corr);
            row.sigipr = static_cast<float>(row.sigipr * corr);
            // LP column = factor the raw intensity was divided by (L·P).
            row.lp = static_cast<float>(1.0 / corr);
        } else {
            row.lp = static_cast<float>(lp_factor(
                request.detector, pm.x_mm, pm.y_mm,
                request.integration.polarization_fraction
            ));
        }

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
            path, crystal, result.mosflm_reflections, wavelength_A,
            result.batch_scans
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
