#include "mosflm/diffraction/StandardProfile.hpp"

#include "mosflm/diffraction/LpCorrection.hpp"
#include "mosflm/geometry/LabFrame.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace mosflm::diffraction {
namespace {

inline std::int32_t pixel_at(
    const image::DiffractionImage& image,
    int col,
    int row
) {
    const int w = image.metadata.width;
    const int h = image.metadata.height;
    if (col < 0 || row < 0 || col >= w || row >= h) {
        return 0;
    }
    return image.pixels[static_cast<std::size_t>(row) * static_cast<std::size_t>(w) +
                        static_cast<std::size_t>(col)];
}

void to_pixel(
    const DetectorGeometry& detector,
    const IntegrationSettings& settings,
    double x_mm,
    double y_mm,
    int& col,
    int& row
) {
    if (settings.mosflm_detector_axes) {
        double xdet = 0.0;
        double ydet = 0.0;
        mosflm::cpp_mm_to_mosflm_xdet_yscale(
            x_mm,
            y_mm,
            settings.image_beam_x_mm,
            settings.image_beam_y_mm,
            detector.pixel_size_x_mm,
            detector.pixel_size_y_mm,
            detector.yscale,
            xdet,
            ydet
        );
        col = static_cast<int>(std::lround(xdet));
        row = static_cast<int>(std::lround(ydet));
        return;
    }
    SummationIntegrator::mm_to_image_pixel(detector, x_mm, y_mm, col, row);
}

int derive_peak_half(
    const DetectorGeometry& detector,
    const MeasurementBoxSettings& box
) {
    if (box.peak_half_px >= 0) {
        return box.peak_half_px;
    }
    const double px = std::max(detector.pixel_size_x_mm, 1e-9);
    const double py = std::max(detector.pixel_size_y_mm, 1e-9);
    const int hx = static_cast<int>(std::lround(0.5 * box.separation_x_mm / px));
    const int hy = static_cast<int>(std::lround(0.5 * box.separation_y_mm / py));
    return std::max(2, std::min(hx, hy));
}

std::vector<double> make_gaussian_profile(int half, double sigma_px) {
    const int n = 2 * half + 1;
    std::vector<double> p(static_cast<std::size_t>(n * n), 0.0);
    const double inv = 1.0 / (2.0 * sigma_px * sigma_px);
    double sum = 0.0;
    for (int dr = -half; dr <= half; ++dr) {
        for (int dc = -half; dc <= half; ++dc) {
            const double w =
                std::exp(-static_cast<double>(dc * dc + dr * dr) * inv);
            p[static_cast<std::size_t>((dr + half) * n + (dc + half))] = w;
            sum += w;
        }
    }
    if (sum > 0.0) {
        for (double& v : p) {
            v /= sum;
        }
    }
    return p;
}

void normalize_profile(std::vector<double>& profile) {
    // Learned averages of signed pixels can dip just below 0 in the tails.
    double sum = 0.0;
    for (double& v : profile) {
        v = std::max(v, 0.0);
        sum += v;
    }
    if (sum > 0.0) {
        for (double& v : profile) {
            v /= sum;
        }
    }
}

struct ProfileMoments {
    double pr1{0.0};
    double pr2{0.0};
    double pr3{0.0};
    double pr4{0.0};
};

ProfileMoments moments_of(
    const std::vector<double>& profile,
    int peak_half,
    int side
) {
    ProfileMoments m;
    for (int dr = -peak_half; dr <= peak_half; ++dr) {
        for (int dc = -peak_half; dc <= peak_half; ++dc) {
            const double P = profile[static_cast<std::size_t>(
                (dr + peak_half) * side + (dc + peak_half)
            )];
            const double p = static_cast<double>(dc);
            const double q = static_cast<double>(dr);
            m.pr1 += p * P;
            m.pr2 += q * P;
            m.pr3 += P;
            m.pr4 += P * P;
        }
    }
    return m;
}

int region_index(
    int col,
    int row,
    int width,
    int height,
    int n_bins_x,
    int n_bins_y
) {
    const int nx = std::max(1, n_bins_x);
    const int ny = std::max(1, n_bins_y);
    const int w = std::max(width, 1);
    const int h = std::max(height, 1);
    int ix = (col * nx) / w;
    int iy = (row * ny) / h;
    if (ix < 0) {
        ix = 0;
    }
    if (iy < 0) {
        iy = 0;
    }
    if (ix >= nx) {
        ix = nx - 1;
    }
    if (iy >= ny) {
        iy = ny - 1;
    }
    return iy * nx + ix;
}

}  // namespace

std::vector<ProfileObservation> StandardProfileFitter::fit(
    const IntegrationRequest& request,
    const StandardProfileSettings& settings
) const {
    IntegrationSettings iset_no_lp = request.settings;
    iset_no_lp.apply_lp = false;
    const IntegrationRequest ireq_meas{
        request.image, request.detector, request.reflections, iset_no_lp
    };
    const auto meas =
        MeasIntegrator{}.integrate(ireq_meas, settings.box);

    const auto& image = request.image;
    const auto& detector = request.detector;
    const auto& predictions = request.reflections;
    const auto& iset = request.settings;
    const int n = static_cast<int>(predictions.size());
    std::vector<ProfileObservation> out;
    out.reserve(static_cast<std::size_t>(n));
    if (n == 0) {
        return out;
    }

    const int peak_half = derive_peak_half(detector, settings.box);
    const int rim = std::max(0, settings.box.bg_rim_px);
    const int outer = peak_half + rim;
    const int side = 2 * peak_half + 1;
    const int np = side * side;
    const int n_bins_x = std::max(1, settings.n_bins_x);
    const int n_bins_y = std::max(1, settings.n_bins_y);
    const int n_regions = n_bins_x * n_bins_y;

    std::vector<int> col(static_cast<std::size_t>(n));
    std::vector<int> row(static_cast<std::size_t>(n));
    std::vector<int> region(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        to_pixel(
            detector, iset, predictions[static_cast<std::size_t>(i)].x_mm,
            predictions[static_cast<std::size_t>(i)].y_mm,
            col[static_cast<std::size_t>(i)], row[static_cast<std::size_t>(i)]
        );
        if (settings.box.refine_centroid) {
            refine_centroid_com(
                image, col[static_cast<std::size_t>(i)],
                row[static_cast<std::size_t>(i)], settings.box.centroid_half_px,
                settings.box.centroid_max_shift_px
            );
        }
        region[static_cast<std::size_t>(i)] = region_index(
            col[static_cast<std::size_t>(i)],
            row[static_cast<std::size_t>(i)],
            image.metadata.width,
            image.metadata.height,
            n_bins_x,
            n_bins_y
        );
    }

    const double sep_x_px =
        settings.box.separation_x_mm / std::max(detector.pixel_size_x_mm, 1e-9);
    const double sep_y_px =
        settings.box.separation_y_mm / std::max(detector.pixel_size_y_mm, 1e-9);
    const double neigh_x = std::max(sep_x_px, static_cast<double>(2 * outer + 1));
    const double neigh_y = std::max(sep_y_px, static_cast<double>(2 * outer + 1));
    std::vector<std::vector<int>> neighbours(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double dx = static_cast<double>(
                col[static_cast<std::size_t>(i)] - col[static_cast<std::size_t>(j)]
            );
            const double dy = static_cast<double>(
                row[static_cast<std::size_t>(i)] - row[static_cast<std::size_t>(j)]
            );
            if (std::fabs(dx) < neigh_x && std::fabs(dy) < neigh_y) {
                neighbours[static_cast<std::size_t>(i)].push_back(j);
                neighbours[static_cast<std::size_t>(j)].push_back(i);
            }
        }
    }

    // Accumulate per-region + global learned profiles.
    std::vector<std::vector<double>> region_acc(
        static_cast<std::size_t>(n_regions),
        std::vector<double>(static_cast<std::size_t>(np), 0.0)
    );
    std::vector<int> region_nlearn(static_cast<std::size_t>(n_regions), 0);
    std::vector<double> global_acc(static_cast<std::size_t>(np), 0.0);
    int global_nlearn = 0;

    for (int i = 0; i < n; ++i) {
        const auto& m = meas[static_cast<std::size_t>(i)];
        // Overloaded spots have clipped peaks: never learn shapes from them.
        if (!m.ok || m.sigma <= 0.0 || m.n_overload > 0) {
            continue;
        }
        if (m.intensity / m.sigma < settings.learn_min_isigi) {
            continue;
        }
        const int c0 = col[static_cast<std::size_t>(i)];
        const int r0 = row[static_cast<std::size_t>(i)];
        const int reg = region[static_cast<std::size_t>(i)];
        std::vector<int> nc, nr;
        for (int j : neighbours[static_cast<std::size_t>(i)]) {
            nc.push_back(col[static_cast<std::size_t>(j)]);
            nr.push_back(row[static_cast<std::size_t>(j)]);
        }
        const auto* tmask =
            settings.box.trusted_mask.empty() ? nullptr
                                              : &settings.box.trusted_mask;
        const BackgroundPlane bg = fit_background_plane(
            image, c0, r0, peak_half, outer, nc, nr,
            settings.box.bg_reject_nsigma, tmask, settings.box.trusted_width
        );
        if (!bg.ok) {
            continue;
        }
        // Signed bg-subtracted values: clipping at 0 turns background noise
        // into a positive pedestal over the whole box, and a profile with a
        // pedestal over-reads sharp spots by ~1/(1 − pedestal fraction).
        // Each learner is normalised to unit volume (not its noisy maximum).
        double volume = 0.0;
        std::vector<double> local(static_cast<std::size_t>(np), 0.0);
        for (int dr = -peak_half; dr <= peak_half; ++dr) {
            for (int dc = -peak_half; dc <= peak_half; ++dc) {
                const int c = c0 + dc;
                const int r = r0 + dr;
                const double od = static_cast<double>(pixel_at(image, c, r));
                const double b = bg.at(
                    static_cast<double>(dc), static_cast<double>(dr)
                );
                const double v = od - b;
                const std::size_t idx = static_cast<std::size_t>(
                    (dr + peak_half) * side + (dc + peak_half)
                );
                local[idx] = v;
                volume += v;
            }
        }
        if (volume < 1.0) {
            continue;
        }
        for (int k = 0; k < np; ++k) {
            const double contrib =
                local[static_cast<std::size_t>(k)] / volume;
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
        global_profile =
            make_gaussian_profile(peak_half, settings.fallback_sigma_px);
    }

    std::vector<std::vector<double>> profiles(
        static_cast<std::size_t>(n_regions)
    );
    std::vector<ProfileMoments> moms(static_cast<std::size_t>(n_regions));
    for (int r = 0; r < n_regions; ++r) {
        if (region_nlearn[static_cast<std::size_t>(r)] >= settings.min_learners) {
            profiles[static_cast<std::size_t>(r)] =
                region_acc[static_cast<std::size_t>(r)];
            normalize_profile(profiles[static_cast<std::size_t>(r)]);
        } else {
            profiles[static_cast<std::size_t>(r)] = global_profile;
        }
        moms[static_cast<std::size_t>(r)] =
            moments_of(profiles[static_cast<std::size_t>(r)], peak_half, side);
    }

    for (int i = 0; i < n; ++i) {
        ProfileObservation obs;
        static_cast<MeasObservation&>(obs) = meas[static_cast<std::size_t>(i)];
        obs.summation_intensity = obs.intensity;
        obs.summation_sigma = obs.sigma;
        obs.profile_region = region[static_cast<std::size_t>(i)];

        if (!obs.ok) {
            out.push_back(obs);
            continue;
        }

        const int c0 = col[static_cast<std::size_t>(i)];
        const int r0 = row[static_cast<std::size_t>(i)];
        const int reg = region[static_cast<std::size_t>(i)];
        const auto& profile = profiles[static_cast<std::size_t>(reg)];
        const auto& pr = moms[static_cast<std::size_t>(reg)];

        std::vector<int> nc, nr;
        for (int j : neighbours[static_cast<std::size_t>(i)]) {
            nc.push_back(col[static_cast<std::size_t>(j)]);
            nr.push_back(row[static_cast<std::size_t>(j)]);
        }
        const auto* tmask2 =
            settings.box.trusted_mask.empty() ? nullptr
                                              : &settings.box.trusted_mask;
        const BackgroundPlane bg = fit_background_plane(
            image, c0, r0, peak_half, outer, nc, nr,
            settings.box.bg_reject_nsigma, tmask2, settings.box.trusted_width
        );
        obs.bg_a = bg.a;
        obs.bg_b = bg.b;
        obs.bg_c = bg.c;

        if (!bg.ok || pr.pr4 < 1e-18) {
            obs.profile_ok = false;
            out.push_back(obs);
            continue;
        }

        // Least-squares profile fit over valid pixels only: overloaded
        // (−2) pixels are skipped, so saturated spots are recovered from
        // their unsaturated wings (MOSFLM profile-fitted overloads).
        double numer = 0.0;
        double pr4_valid = 0.0;
        double p_valid = 0.0;   // Σ P
        double p2b_valid = 0.0; // Σ P² b
        double p3_valid = 0.0;  // Σ P³
        double bg_under = 0.0;
        int n_peak = 0;
        for (int dr = -peak_half; dr <= peak_half; ++dr) {
            for (int dc = -peak_half; dc <= peak_half; ++dc) {
                const int c = c0 + dc;
                const int rr = r0 + dr;
                const std::int32_t raw = pixel_at(image, c, rr);
                const double b = bg.a * static_cast<double>(dc) +
                                 bg.b * static_cast<double>(dr) + bg.c;
                // Background plane under the peak (for Poisson σ model).
                bg_under += b;
                ++n_peak;
                if (raw < 0) {
                    continue;
                }
                const double P = profile[static_cast<std::size_t>(
                    (dr + peak_half) * side + (dc + peak_half)
                )];
                numer += P * (static_cast<double>(raw) - b);
                pr4_valid += P * P;
                p_valid += P;
                p2b_valid += P * P * std::max(b, 0.0);
                p3_valid += P * P * P;
            }
        }
        if (pr4_valid < 1e-18) {
            obs.profile_ok = false;
            out.push_back(obs);
            continue;
        }
        const double iprf = numer * pr.pr3 / pr4_valid;
        obs.profile_intensity = iprf;
        // Poisson variance of the weighted LS fit (model pixel b + k·P):
        //   var(k) = Σ P² (b + k P) / (Σ P²)²,   I = k · ΣP_all
        // plus the error of the background estimate (mean over n_bg rim
        // pixels) propagated through the fit: δk = −δb · ΣP / ΣP².
        const double bg_mean =
            (n_peak > 0) ? (bg_under / static_cast<double>(n_peak)) : 0.0;
        const int n_bg = std::max(1, obs.npix_bg);
        // Signed k keeps E[var] unbiased for weak spots; max(k, 0) overstates
        // σ when the true I is ~0.
        const double k_fit = numer / pr4_valid;
        const double var_k_peak =
            std::max(p2b_valid + k_fit * p3_valid, 0.0) /
            (pr4_valid * pr4_valid);
        const double dk_db = p_valid / pr4_valid;
        const double var_k_bg =
            dk_db * dk_db * std::fabs(bg_mean) / static_cast<double>(n_bg);
        double var = pr.pr3 * pr.pr3 * (var_k_peak + var_k_bg);
        // No unit floor: a sparse Eiger slice can have var < 1, and partial
        // slices are summed (≈7 per reflection), so a floor of 1 per slice
        // inflates σ of weak reflections. Keep σ > 0 only.
        var = std::max(var, 1e-6);
        obs.profile_sigma = std::sqrt(var);
        obs.pkratio =
            std::fabs(obs.summation_intensity - iprf) /
            std::max(std::fabs(iprf), 1.0);
        obs.profile_ok = true;

        if (settings.prefer_profile &&
            obs.pkratio <= settings.max_pkratio_for_profile) {
            obs.intensity = obs.profile_intensity;
            obs.sigma = obs.profile_sigma;
        }
        if (iset.apply_lp) {
            const double rho = iset.polarization_fraction;
            apply_lp_inplace(
                detector, obs.x_mm, obs.y_mm, obs.intensity, obs.sigma, rho
            );
            apply_lp_inplace(
                detector, obs.x_mm, obs.y_mm, obs.profile_intensity,
                obs.profile_sigma, rho
            );
            apply_lp_inplace(
                detector, obs.x_mm, obs.y_mm, obs.summation_intensity,
                obs.summation_sigma, rho
            );
        }
        out.push_back(obs);
    }
    return out;
}

}  // namespace mosflm::diffraction
