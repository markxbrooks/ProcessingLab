#include "mosflm/diffraction/Profile3d.hpp"

#include "mosflm/diffraction/LpCorrection.hpp"
#include "mosflm/diffraction/TrustedRegion.hpp"
#include "mosflm/geometry/LabFrame.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace mosflm::diffraction {
namespace {

constexpr double kPi = 3.14159265358979323846;

bool voxel_trusted(
    const Profile3dSettings& settings,
    const Shoebox& box,
    int dc,
    int dr
) {
    if (settings.box.trusted_mask.empty() || settings.box.trusted_width <= 0) {
        return true;
    }
    return is_trusted(
        settings.box.trusted_mask,
        settings.box.trusted_width,
        box.col0 + dc,
        box.row0 + dr
    );
}

int region_of(
    int col,
    int row,
    int width,
    int height,
    int n_bins_x,
    int n_bins_y
) {
    const int nx = std::max(1, n_bins_x);
    const int ny = std::max(1, n_bins_y);
    const int w = std::max(1, width);
    const int h = std::max(1, height);
    int ix = (col * nx) / w;
    int iy = (row * ny) / h;
    ix = std::clamp(ix, 0, nx - 1);
    iy = std::clamp(iy, 0, ny - 1);
    return iy * nx + ix;
}

void normalize_profile(std::vector<double>& p) {
    double s = 0.0;
    for (double v : p) {
        s += v;
    }
    if (!(s > 1e-18)) {
        return;
    }
    const double inv = 1.0 / s;
    for (double& v : p) {
        v *= inv;
    }
}

std::vector<double> make_gaussian_3d(
    int peak_half,
    int nphi,
    double sig_xy,
    double sig_phi
) {
    const int side = 2 * peak_half + 1;
    std::vector<double> p(
        static_cast<std::size_t>(nphi) * static_cast<std::size_t>(side) *
        static_cast<std::size_t>(side),
        0.0
    );
    const double sxy = std::max(0.25, sig_xy);
    const double sph = std::max(0.25, sig_phi);
    const double inv_xy = 1.0 / (2.0 * sxy * sxy);
    const double inv_ph = 1.0 / (2.0 * sph * sph);
    const int cphi = nphi / 2;
    for (int iphi = 0; iphi < nphi; ++iphi) {
        const double dp = static_cast<double>(iphi - cphi);
        const double gphi = std::exp(-dp * dp * inv_ph);
        for (int dr = -peak_half; dr <= peak_half; ++dr) {
            for (int dc = -peak_half; dc <= peak_half; ++dc) {
                const double r2 =
                    static_cast<double>(dc * dc + dr * dr);
                const std::size_t idx =
                    static_cast<std::size_t>(iphi) *
                        static_cast<std::size_t>(side * side) +
                    static_cast<std::size_t>(dr + peak_half) *
                        static_cast<std::size_t>(side) +
                    static_cast<std::size_t>(dc + peak_half);
                p[idx] = std::exp(-r2 * inv_xy) * gphi;
            }
        }
    }
    normalize_profile(p);
    return p;
}

std::size_t peak_index(int iphi, int dr, int dc, int peak_half, int nphi) {
    const int side = 2 * peak_half + 1;
    (void)nphi;
    return static_cast<std::size_t>(iphi) * static_cast<std::size_t>(side * side) +
           static_cast<std::size_t>(dr + peak_half) * static_cast<std::size_t>(side) +
           static_cast<std::size_t>(dc + peak_half);
}

Background3d fit_background_3d(
    const Shoebox& box,
    const Profile3dSettings& settings
) {
    Background3d out;
    const int peak_half = box.peak_half;
    const int outer = box.outer_half;
    const int nphi = box.nphi;
    const int cphi = nphi / 2;
    const bool use_phi = settings.use_phi_bg;
    const int npar = use_phi ? 4 : 3;
    const double nsigma = settings.box.bg_reject_nsigma;

    struct Sample {
        double dc{};
        double dr{};
        double dphi{};
        double od{};
        bool keep{true};
    };
    std::vector<Sample> samples;
    samples.reserve(static_cast<std::size_t>(nphi) * 64);

    for (int iphi = 0; iphi < nphi; ++iphi) {
        const double dphi = static_cast<double>(iphi - cphi);
        for (int dr = -outer; dr <= outer; ++dr) {
            for (int dc = -outer; dc <= outer; ++dc) {
                const int adc = std::abs(dc);
                const int adr = std::abs(dr);
                if (adc <= peak_half && adr <= peak_half) {
                    continue;
                }
                if (!box.is_valid(iphi, dr, dc)) {
                    continue;
                }
                if (!voxel_trusted(settings, box, dc, dr)) {
                    continue;
                }
                Sample s;
                s.dc = static_cast<double>(dc);
                s.dr = static_cast<double>(dr);
                s.dphi = dphi;
                s.od = static_cast<double>(box.at(iphi, dr, dc));
                samples.push_back(s);
            }
        }
    }
    if (static_cast<int>(samples.size()) < npar + 3) {
        return out;
    }

    auto solve = [&](std::vector<Sample>& pts) -> bool {
        int nkeep = 0;
        for (const auto& s : pts) {
            if (s.keep) {
                ++nkeep;
            }
        }
        if (nkeep < npar + 3) {
            return false;
        }
        Eigen::MatrixXd A(nkeep, npar);
        Eigen::VectorXd b(nkeep);
        int row = 0;
        for (const auto& s : pts) {
            if (!s.keep) {
                continue;
            }
            A(row, 0) = s.dc;
            A(row, 1) = s.dr;
            A(row, 2) = 1.0;
            if (use_phi) {
                A(row, 3) = s.dphi;
            }
            b(row) = s.od;
            ++row;
        }
        const Eigen::VectorXd x =
            A.colPivHouseholderQr().solve(b);
        if (x.size() < npar || !x.allFinite()) {
            return false;
        }
        out.a = x(0);
        out.b = x(1);
        out.c = x(2);
        out.d = use_phi ? x(3) : 0.0;
        out.n_used = nkeep;
        return true;
    };

    if (!solve(samples)) {
        return out;
    }
    if (nsigma > 0.0) {
        for (int pass = 0; pass < 3; ++pass) {
            double mean = 0.0;
            double m2 = 0.0;
            int n = 0;
            for (const auto& s : samples) {
                if (!s.keep) {
                    continue;
                }
                const double r = s.od - out.at(s.dc, s.dr, s.dphi);
                ++n;
                const double d = r - mean;
                mean += d / static_cast<double>(n);
                m2 += d * (r - mean);
            }
            if (n < npar + 3) {
                break;
            }
            const double sd = std::sqrt(std::max(m2 / static_cast<double>(n), 1e-6));
            int dropped = 0;
            for (auto& s : samples) {
                if (!s.keep) {
                    continue;
                }
                // High outliers only: a symmetric clip on sparse Eiger
                // counts rejects every 1-count voxel and zeroes the bg.
                const double model = out.at(s.dc, s.dr, s.dphi);
                const double lim = std::max(
                    model + mean + nsigma * sd,
                    poisson_upper_limit(model, nsigma)
                );
                if (s.od > lim) {
                    s.keep = false;
                    ++dropped;
                }
            }
            if (dropped == 0) {
                break;
            }
            if (!solve(samples)) {
                break;
            }
        }
    }
    out.ok = out.n_used >= npar + 3;
    return out;
}

double phi_window_fraction(int nphi, double sig_phi_frames) {
    const double sph = std::max(0.25, sig_phi_frames);
    const int half = nphi / 2;
    const int span = std::max(20, half + 8);
    double all = 0.0;
    double win = 0.0;
    for (int k = -span; k <= span; ++k) {
        const double w =
            std::exp(-0.5 * static_cast<double>(k * k) / (sph * sph));
        all += w;
        if (k >= -half && k <= half) {
            win += w;
        }
    }
    if (!(all > 0.0)) {
        return 1.0;
    }
    return std::min(1.0, win / all);
}

double fallback_sigma_phi(const Profile3dSettings& settings) {
    if (settings.fallback_sigma_phi_frames > 0.0) {
        const double osc = std::max(1e-6, settings.osc_deg);
        const double mosaic = std::max(0.0, settings.mosaic_spread_deg);
        if (mosaic > 0.0) {
            // σ_M in degrees; at ζ ≈ 1 the frame width is σ_M / osc.
            return std::max(0.4, mosaic / osc);
        }
        return settings.fallback_sigma_phi_frames;
    }
    return 1.0;
}

Eigen::Matrix3d rotation_about_x_deg(double phi_deg) {
    const double r = phi_deg * (kPi / 180.0);
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
    R(1, 1) = std::cos(r);
    R(1, 2) = -std::sin(r);
    R(2, 1) = std::sin(r);
    R(2, 2) = std::cos(r);
    return R;
}

}  // namespace

std::vector<Profile3dObservation> Profile3dFitter::fit(
    const std::vector<Shoebox>& boxes,
    const DetectorGeometry& detector,
    const IntegrationSettings& integration,
    const Profile3dSettings& settings
) const {
    (void)detector;
    (void)integration;
    std::vector<Profile3dObservation> out;
    if (boxes.empty()) {
        return out;
    }

    const int n_bins_x = std::max(1, settings.n_bins_x);
    const int n_bins_y = std::max(1, settings.n_bins_y);
    const int n_regions = n_bins_x * n_bins_y;
    const double sig_phi = fallback_sigma_phi(settings);

    // Pass 1: summation I/σ with a joint 3D background (for learning).
    struct Prep {
        Background3d bg;
        double i_sum{0.0};
        double sig_sum{0.0};
        int npix_peak{0};
        int npix_bg{0};
        double bg_mean{0.0};
        bool ok{false};
    };
    std::vector<Prep> prep(boxes.size());
    int peak_half = 0;
    int nphi = 0;
    for (std::size_t i = 0; i < boxes.size(); ++i) {
        const auto& box = boxes[i];
        Profile3dObservation obs;
        obs.h = box.h;
        obs.k = box.k;
        obs.l = box.l;
        obs.x_mm = box.x_mm;
        obs.y_mm = box.y_mm;
        obs.phi_calc_deg = box.phi_calc_deg;
        obs.peak_frame = box.peak_frame;
        obs.col0 = box.col0;
        obs.row0 = box.row0;
        obs.nphi = box.nphi;
        if (!box.ok) {
            out.push_back(obs);
            continue;
        }
        peak_half = box.peak_half;
        nphi = box.nphi;
        auto& pr = prep[i];
        pr.bg = fit_background_3d(box, settings);
        if (!pr.bg.ok) {
            out.push_back(obs);
            continue;
        }
        const int cphi = box.nphi / 2;
        double i_sum = 0.0;
        double bg_under = 0.0;
        int n_peak = 0;
        int n_bg = pr.bg.n_used;
        for (int iphi = 0; iphi < box.nphi; ++iphi) {
            const double dphi = static_cast<double>(iphi - cphi);
            for (int dr = -box.peak_half; dr <= box.peak_half; ++dr) {
                for (int dc = -box.peak_half; dc <= box.peak_half; ++dc) {
                    if (!box.is_valid(iphi, dr, dc) ||
                        !voxel_trusted(settings, box, dc, dr)) {
                        continue;
                    }
                    const double od = static_cast<double>(box.at(iphi, dr, dc));
                    const double b = pr.bg.at(
                        static_cast<double>(dc), static_cast<double>(dr), dphi
                    );
                    i_sum += od - b;
                    bg_under += b;
                    ++n_peak;
                }
            }
        }
        pr.npix_peak = n_peak;
        pr.npix_bg = std::max(1, n_bg);
        pr.i_sum = i_sum;
        pr.bg_mean =
            (n_peak > 0) ? (bg_under / static_cast<double>(n_peak)) : 0.0;
        // Poisson: all counts in the peak (signal + background under it),
        // plus the background-level error scaled to the peak volume.
        double var = std::max(i_sum + bg_under, 0.0);
        if (n_peak > 0) {
            var += (static_cast<double>(n_peak * n_peak) /
                    static_cast<double>(pr.npix_bg)) *
                   std::fabs(pr.bg_mean);
        }
        pr.sig_sum = std::sqrt(std::max(var, 1.0));
        pr.ok = n_peak > 0;
        obs.ok = pr.ok;
        obs.npix_peak = n_peak;
        obs.npix_bg = pr.npix_bg;
        obs.summation_intensity = i_sum;
        obs.summation_sigma = pr.sig_sum;
        obs.intensity = i_sum;
        obs.sigma = pr.sig_sum;
        out.push_back(obs);
    }

    if (peak_half < 1 || nphi < 1) {
        return out;
    }
    const int side = 2 * peak_half + 1;
    const int np = nphi * side * side;

    std::vector<std::vector<double>> region_acc(
        static_cast<std::size_t>(n_regions),
        std::vector<double>(static_cast<std::size_t>(np), 0.0)
    );
    std::vector<int> region_nlearn(static_cast<std::size_t>(n_regions), 0);
    std::vector<double> global_acc(static_cast<std::size_t>(np), 0.0);
    int global_nlearn = 0;

    for (std::size_t i = 0; i < boxes.size(); ++i) {
        const auto& box = boxes[i];
        const auto& pr = prep[i];
        if (!box.ok || !pr.ok || pr.sig_sum <= 0.0) {
            continue;
        }
        if (pr.i_sum / pr.sig_sum < settings.learn_min_isigi) {
            continue;
        }
        if (box.peak_half != peak_half || box.nphi != nphi) {
            continue;
        }
        const int cphi = nphi / 2;
        std::vector<double> local(static_cast<std::size_t>(np), 0.0);
        double peak_max = 0.0;
        for (int iphi = 0; iphi < nphi; ++iphi) {
            const double dphi = static_cast<double>(iphi - cphi);
            for (int dr = -peak_half; dr <= peak_half; ++dr) {
                for (int dc = -peak_half; dc <= peak_half; ++dc) {
                    if (!box.is_valid(iphi, dr, dc)) {
                        continue;
                    }
                    const double od = static_cast<double>(box.at(iphi, dr, dc));
                    const double b = pr.bg.at(
                        static_cast<double>(dc), static_cast<double>(dr), dphi
                    );
                    const double v = std::max(0.0, od - b);
                    const std::size_t idx =
                        peak_index(iphi, dr, dc, peak_half, nphi);
                    local[idx] = v;
                    peak_max = std::max(peak_max, v);
                }
            }
        }
        if (peak_max < 1.0) {
            continue;
        }
        const int reg = region_of(
            box.col0, box.row0, box.width, box.height, n_bins_x, n_bins_y
        );
        for (int k = 0; k < np; ++k) {
            const double contrib =
                local[static_cast<std::size_t>(k)] / peak_max;
            region_acc[static_cast<std::size_t>(reg)]
                      [static_cast<std::size_t>(k)] += contrib;
            global_acc[static_cast<std::size_t>(k)] += contrib;
        }
        ++region_nlearn[static_cast<std::size_t>(reg)];
        ++global_nlearn;
    }

    std::vector<double> global_profile;
    if (global_nlearn >= settings.min_learners) {
        global_profile = global_acc;
        normalize_profile(global_profile);
    } else {
        global_profile = make_gaussian_3d(
            peak_half, nphi, settings.fallback_sigma_xy_px, sig_phi
        );
    }
    std::vector<std::vector<double>> profiles(
        static_cast<std::size_t>(n_regions)
    );
    for (int r = 0; r < n_regions; ++r) {
        if (region_nlearn[static_cast<std::size_t>(r)] >= settings.min_learners) {
            profiles[static_cast<std::size_t>(r)] =
                region_acc[static_cast<std::size_t>(r)];
            normalize_profile(profiles[static_cast<std::size_t>(r)]);
        } else {
            profiles[static_cast<std::size_t>(r)] = global_profile;
        }
    }

    const double frac_win = phi_window_fraction(nphi, sig_phi);

    // Pass 2: Rossmann 3D fit.
    for (std::size_t i = 0; i < boxes.size(); ++i) {
        auto& obs = out[i];
        const auto& box = boxes[i];
        const auto& pr = prep[i];
        obs.fractioncalc = frac_win;
        if (!obs.ok || !pr.ok || !pr.bg.ok) {
            continue;
        }
        if (box.peak_half != peak_half || box.nphi != nphi) {
            continue;
        }
        const int reg = region_of(
            box.col0, box.row0, box.width, box.height, n_bins_x, n_bins_y
        );
        obs.profile_region = reg;
        const auto& profile = profiles[static_cast<std::size_t>(reg)];
        const int cphi = nphi / 2;

        double pr_dc = 0.0;
        double pr_dr = 0.0;
        double pr_1 = 0.0;
        double pr_p2 = 0.0;
        double pr_dphi = 0.0;
        double sumprod = 0.0;
        double bg_under = 0.0;
        double pr_p2b = 0.0;
        double pr_p3 = 0.0;
        int n_peak = 0;
        for (int iphi = 0; iphi < nphi; ++iphi) {
            const double dphi = static_cast<double>(iphi - cphi);
            for (int dr = -peak_half; dr <= peak_half; ++dr) {
                for (int dc = -peak_half; dc <= peak_half; ++dc) {
                    const std::size_t idx =
                        peak_index(iphi, dr, dc, peak_half, nphi);
                    const double P = profile[idx];
                    if (!box.is_valid(iphi, dr, dc) ||
                        !voxel_trusted(settings, box, dc, dr)) {
                        continue;
                    }
                    const double od = static_cast<double>(box.at(iphi, dr, dc));
                    sumprod += P * od;
                    pr_dc += P * static_cast<double>(dc);
                    pr_dr += P * static_cast<double>(dr);
                    pr_1 += P;
                    pr_p2 += P * P;
                    pr_dphi += P * dphi;
                    const double b = pr.bg.at(
                        static_cast<double>(dc), static_cast<double>(dr), dphi
                    );
                    bg_under += b;
                    pr_p2b += P * P * std::max(b, 0.0);
                    pr_p3 += P * P * P;
                    ++n_peak;
                }
            }
        }
        if (!(pr_p2 > 1e-18) || !(pr_1 > 1e-18)) {
            continue;
        }
        const double numer =
            sumprod - pr.bg.a * pr_dc - pr.bg.b * pr_dr - pr.bg.c * pr_1 -
            pr.bg.d * pr_dphi;
        const double iprf = numer * pr_1 / pr_p2;
        obs.profile_intensity = iprf;
        const double bg_mean =
            (n_peak > 0) ? (bg_under / static_cast<double>(n_peak)) : 0.0;
        const int n_bg = std::max(1, obs.npix_bg);
        // I = k·ΣP(od − b) with k = ΣP/ΣP²: each voxel holds Poisson counts
        // m = b + I·P/ΣP, and the fitted background level adds
        // (k·ΣP)²·bg/n_bg (n_peak²·bg/n_bg for a flat profile). The signed
        // I keeps E[var] unbiased for weak spots (as raw counts do for
        // summation); max(I, 0) would overstate σ when the true I is ~0.
        const double k = pr_1 / pr_p2;
        double var = k * k * (pr_p2b + iprf / pr_1 * pr_p3);
        if (n_peak > 0) {
            var += k * k * pr_1 * pr_1 / static_cast<double>(n_bg) *
                   std::fabs(bg_mean);
        }
        obs.profile_sigma = std::sqrt(std::max(var, 1.0));
        obs.pkratio =
            std::fabs(obs.summation_intensity - iprf) /
            std::max(std::fabs(iprf), 1.0);
        obs.profile_ok = true;
        if (settings.prefer_profile &&
            obs.pkratio <= settings.max_pkratio_for_profile) {
            obs.intensity = obs.profile_intensity;
            obs.sigma = obs.profile_sigma;
        }
    }
    return out;
}

crystal::MosflmReflection profile3d_to_mosflm(
    const Profile3dObservation& obs,
    const DetectorGeometry& detector,
    const OrientationMatrix& orientation,
    const crystal::CrystalMetadata& crystal,
    const IntegrationSettings& integration,
    double air_absorption_mm_inv,
    double osc_deg
) {
    crystal::MosflmReflection row;
    row.h = obs.h;
    row.k = obs.k;
    row.l = obs.l;
    row.h0 = obs.h;
    row.k0 = obs.k;
    row.l0 = obs.l;
    row.batch = std::max(1, obs.peak_frame);
    // MOSFLM convention: I = summation, IPR = profile (as ProcessPipeline).
    row.intensity = static_cast<float>(obs.summation_intensity);
    row.sigma = static_cast<float>(obs.summation_sigma);
    row.ipr = static_cast<float>(
        obs.profile_ok ? obs.profile_intensity : obs.summation_intensity
    );
    row.sigipr = static_cast<float>(
        obs.profile_ok ? obs.profile_sigma : obs.summation_sigma
    );
    row.fractioncalc = static_cast<float>(obs.fractioncalc);
    row.mpart = obs.fractioncalc < 0.99 ? 1 : 0;
    row.rot = static_cast<float>(obs.phi_calc_deg);
    row.phi_obs = static_cast<float>(obs.phi_calc_deg);
    const double width =
        osc_deg > 0.0
            ? osc_deg * static_cast<double>(std::max(1, obs.nphi))
            : static_cast<double>(std::max(1, obs.nphi));
    row.width = static_cast<float>(width);

    double xdet = 0.0;
    double ydet = 0.0;
    const double beam_x = integration.mosflm_detector_axes
                              ? integration.image_beam_x_mm
                              : detector.beam_x_mm;
    const double beam_y = integration.mosflm_detector_axes
                              ? integration.image_beam_y_mm
                              : detector.beam_y_mm;
    cpp_mm_to_mosflm_xdet_yscale(
        obs.x_mm,
        obs.y_mm,
        beam_x,
        beam_y,
        detector.pixel_size_x_mm,
        detector.pixel_size_y_mm,
        detector.yscale,
        xdet,
        ydet
    );
    xdet = static_cast<double>(obs.col0);
    ydet = static_cast<double>(obs.row0);
    row.xdet = static_cast<float>(xdet);
    row.ydet = static_cast<float>(ydet);
    row.xdet_obs = row.xdet;
    row.ydet_obs = row.ydet;

    if (integration.rotation_lp) {
        const double lam = detector.wavelength_A;
        const Eigen::Vector3d s0(0.0, 0.0, lam > 0.0 ? 1.0 / lam : 1.0);
        const double mid = std::isfinite(obs.phi_calc_deg) ? obs.phi_calc_deg : 0.0;
        const Eigen::Vector3d s1 =
            s0 + rotation_about_x_deg(mid) *
                     orientation.miller_to_lab(obs.h, obs.k, obs.l);
        const double corr = rotation_lp_correction(
            s0,
            s1,
            Eigen::Vector3d::UnitX(),
            Eigen::Vector3d::UnitY(),
            integration.polarization_fraction
        );
        row.intensity = static_cast<float>(row.intensity * corr);
        row.sigma = static_cast<float>(row.sigma * corr);
        row.ipr = static_cast<float>(row.ipr * corr);
        row.sigipr = static_cast<float>(row.sigipr * corr);
        row.lp = static_cast<float>(1.0 / corr);
    } else {
        row.lp = static_cast<float>(lp_factor(
            detector, obs.x_mm, obs.y_mm, integration.polarization_fraction
        ));
    }

    if (air_absorption_mm_inv > 0.0 && detector.distance_mm > 0.0) {
        const double dx = obs.x_mm - detector.beam_x_mm;
        const double dy = obs.y_mm - detector.beam_y_mm;
        const double rad = std::hypot(dx, dy);
        const double path = std::hypot(detector.distance_mm, rad);
        const double extra = std::max(0.0, path - detector.distance_mm);
        const double corr = std::exp(air_absorption_mm_inv * extra);
        row.intensity = static_cast<float>(row.intensity * corr);
        row.sigma = static_cast<float>(row.sigma * corr);
        row.ipr = static_cast<float>(row.ipr * corr);
        row.sigipr = static_cast<float>(row.sigipr * corr);
    }

    if (crystal.spacegroup != nullptr) {
        crystal::MtzAdapter::reduce_to_asu(
            crystal.spacegroup, row.h, row.k, row.l, row.misym
        );
    }
    return row;
}

}  // namespace mosflm::diffraction
