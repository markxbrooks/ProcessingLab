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

/// Predicted spot centre in fractional image pixels (column, row).
void to_pixel_f(
    const DetectorGeometry& detector,
    const IntegrationSettings& settings,
    double x_mm,
    double y_mm,
    double& col,
    double& row
) {
    if (settings.mosflm_detector_axes) {
        mosflm::cpp_mm_to_mosflm_xdet_yscale(
            x_mm, y_mm, settings.image_beam_x_mm, settings.image_beam_y_mm,
            detector.pixel_size_x_mm, detector.pixel_size_y_mm, detector.yscale,
            col, row
        );
        return;
    }
    double xf = 0.0;
    double yf = 0.0;
    detector.mm_to_pixel(x_mm, y_mm, xf, yf);
    col = xf + 0.5 * detector.nfast;
    row = yf + 0.5 * detector.nslow;
}

/// Oversampled profile: value of a pixel (fraction of the spot volume) as a
/// function of that pixel centre's offset (px) from the predicted centre,
/// on a grid of ``s`` cells per pixel spanning ±(half + 0.5) px.
struct FineGrid {
    int half{0};
    int s{1};
    int side{1};  ///< Cells per row: (2·half + 1)·s.

    FineGrid(int h, int over) : half(h), s(std::max(over, 1)), side((2 * h + 1) * std::max(over, 1)) {}

    [[nodiscard]] std::size_t size() const {
        return static_cast<std::size_t>(side) * static_cast<std::size_t>(side);
    }
    /// Continuous cell coordinate of offset ``x`` (px).
    [[nodiscard]] double cell(double x) const {
        return (x + half + 0.5) * s - 0.5;
    }
    /// Bilinear sample of the profile (sum/w) at offset (x, y), skipping
    /// empty cells; negatives clipped (learned tails dip below 0).
    [[nodiscard]] double sample(const double* sum, const double* w, double x, double y) const {
        const double u = cell(x);
        const double t = cell(y);
        const int i0 = static_cast<int>(std::floor(u));
        const int j0 = static_cast<int>(std::floor(t));
        const double fu = u - i0;
        const double ft = t - j0;
        double acc = 0.0;
        double wsum = 0.0;
        for (int dj = 0; dj <= 1; ++dj) {
            for (int di = 0; di <= 1; ++di) {
                const int i = std::clamp(i0 + di, 0, side - 1);
                const int j = std::clamp(j0 + dj, 0, side - 1);
                const auto k = static_cast<std::size_t>(j * side + i);
                if (!(w[k] > 0.0)) {
                    continue;
                }
                const double wt = (di ? fu : 1.0 - fu) * (dj ? ft : 1.0 - ft);
                acc += wt * sum[k] / w[k];
                wsum += wt;
            }
        }
        return wsum > 0.0 ? std::max(acc / wsum, 0.0) : 0.0;
    }
};

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

void ProfileAccumulator::add(const ProfileAccumulator& other) {
    if (other.empty()) {
        return;
    }
    if (empty() && region_sum.empty()) {
        *this = other;
        return;
    }
    if (other.half != half || other.oversample != oversample ||
        other.n_bins_x != n_bins_x || other.n_bins_y != n_bins_y ||
        other.region_sum.size() != region_sum.size() ||
        other.global_sum.size() != global_sum.size()) {
        return;
    }
    for (std::size_t k = 0; k < region_sum.size(); ++k) {
        region_sum[k] += other.region_sum[k];
        region_w[k] += other.region_w[k];
    }
    for (std::size_t r = 0; r < region_n.size(); ++r) {
        region_n[r] += other.region_n[r];
    }
    for (std::size_t k = 0; k < global_sum.size(); ++k) {
        global_sum[k] += other.global_sum[k];
        global_w[k] += other.global_w[k];
    }
    global_n += other.global_n;
}

std::vector<ProfileObservation> StandardProfileFitter::fit(
    const IntegrationRequest& request,
    const StandardProfileSettings& settings
) const {
    return fit(request, settings, nullptr);
}

std::vector<ProfileObservation> StandardProfileFitter::fit(
    const IntegrationRequest& request,
    const StandardProfileSettings& settings,
    ProfileAccumulator* learned
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

    // XDS-style: every spot is learned and fitted at its predicted
    // sub-pixel position (no per-spot re-centring), so positional scatter
    // blurs learners and fitted spots alike; the oversampled profile keeps
    // sharp spots from being fitted with a sub-pixel-averaged (broader)
    // shape. Integer centring read thaumatin's ~1 px Eiger spots 15% high
    // at low resolution and its weak high-resolution spots ~45% low.
    const int S = std::max(1, settings.oversample);
    const FineGrid grid(peak_half, S);
    const std::size_t nf = grid.size();

    std::vector<int> col(static_cast<std::size_t>(n));
    std::vector<int> row(static_cast<std::size_t>(n));
    std::vector<double> fx(static_cast<std::size_t>(n));
    std::vector<double> fy(static_cast<std::size_t>(n));
    std::vector<int> region(static_cast<std::size_t>(n));
    // Rings of 2θ (8° wide) after the tiles: spot size varies mainly with
    // radius, so a tile short of learners falls back to its ring before the
    // global profile (dominated by sharp low-angle spots).
    constexpr int kRings = 6;
    std::vector<int> ring(static_cast<std::size_t>(n));
    const int n_slots = n_regions + kRings;
    for (int i = 0; i < n; ++i) {
        const auto ui = static_cast<std::size_t>(i);
        double xc = 0.0;
        double yc = 0.0;
        to_pixel_f(detector, iset, predictions[ui].x_mm, predictions[ui].y_mm, xc, yc);
        col[ui] = static_cast<int>(std::lround(xc));
        row[ui] = static_cast<int>(std::lround(yc));
        fx[ui] = xc - col[ui];
        fy[ui] = yc - row[ui];
        region[ui] = region_index(
            col[ui], row[ui], image.metadata.width, image.metadata.height,
            n_bins_x, n_bins_y
        );
        const double tth = detector.two_theta_deg(predictions[ui].x_mm, predictions[ui].y_mm);
        ring[ui] = n_regions + std::clamp(static_cast<int>(tth / 8.0), 0, kRings - 1);
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
    const auto* tmask =
        settings.box.trusted_mask.empty() ? nullptr : &settings.box.trusted_mask;
    const auto background = [&](int i) {
        std::vector<int> nc;
        std::vector<int> nr;
        for (int j : neighbours[static_cast<std::size_t>(i)]) {
            nc.push_back(col[static_cast<std::size_t>(j)]);
            nr.push_back(row[static_cast<std::size_t>(j)]);
        }
        return fit_background_plane(
            image, col[static_cast<std::size_t>(i)], row[static_cast<std::size_t>(i)],
            peak_half, outer, nc, nr, settings.box.bg_reject_nsigma, tmask,
            settings.box.trusted_width
        );
    };

    // Accumulate per-region + global oversampled profiles.
    std::vector<double> region_sum(static_cast<std::size_t>(n_slots) * nf, 0.0);
    std::vector<double> region_w(static_cast<std::size_t>(n_slots) * nf, 0.0);
    std::vector<int> region_nlearn(static_cast<std::size_t>(n_slots), 0);
    std::vector<double> global_sum(nf, 0.0);
    std::vector<double> global_w(nf, 0.0);
    int global_nlearn = 0;

    std::vector<double> local(static_cast<std::size_t>(np), 0.0);
    for (int i = 0; i < n; ++i) {
        const auto ui = static_cast<std::size_t>(i);
        const auto& m = meas[ui];
        // Overloaded spots have clipped peaks: never learn shapes from them.
        if (!m.ok || m.sigma <= 0.0 || m.n_overload > 0) {
            continue;
        }
        if (m.intensity / m.sigma < settings.learn_min_isigi) {
            continue;
        }
        const BackgroundPlane bg = background(i);
        if (!bg.ok) {
            continue;
        }
        // Signed bg-subtracted values: clipping at 0 turns background noise
        // into a positive pedestal over the whole box, and a profile with a
        // pedestal over-reads sharp spots by ~1/(1 − pedestal fraction).
        // Each learner is normalised to unit volume (not its noisy maximum).
        double volume = 0.0;
        for (int dr = -peak_half; dr <= peak_half; ++dr) {
            for (int dc = -peak_half; dc <= peak_half; ++dc) {
                const double od = static_cast<double>(pixel_at(image, col[ui] + dc, row[ui] + dr));
                const double v = od - bg.at(static_cast<double>(dc), static_cast<double>(dr));
                local[static_cast<std::size_t>((dr + peak_half) * side + (dc + peak_half))] = v;
                volume += v;
            }
        }
        if (volume < 1.0) {
            continue;
        }
        const std::size_t base = static_cast<std::size_t>(region[ui]) * nf;
        const std::size_t rbase = static_cast<std::size_t>(ring[ui]) * nf;
        for (int dr = -peak_half; dr <= peak_half; ++dr) {
            for (int dc = -peak_half; dc <= peak_half; ++dc) {
                const double v =
                    local[static_cast<std::size_t>((dr + peak_half) * side + (dc + peak_half))] /
                    volume;
                const double x = dc - fx[ui];
                const double y = dr - fy[ui];
                // Region and global grids share the splat weights.
                const double u = grid.cell(x);
                const double t = grid.cell(y);
                const int i0 = static_cast<int>(std::floor(u));
                const int j0 = static_cast<int>(std::floor(t));
                for (int dj = 0; dj <= 1; ++dj) {
                    for (int di = 0; di <= 1; ++di) {
                        const int ci = i0 + di;
                        const int cj = j0 + dj;
                        if (ci < 0 || cj < 0 || ci >= grid.side || cj >= grid.side) {
                            continue;
                        }
                        const double wt = (di ? u - i0 : 1.0 - (u - i0)) *
                                          (dj ? t - j0 : 1.0 - (t - j0));
                        const auto k = static_cast<std::size_t>(cj * grid.side + ci);
                        region_sum[base + k] += wt * v;
                        region_w[base + k] += wt;
                        region_sum[rbase + k] += wt * v;
                        region_w[rbase + k] += wt;
                        global_sum[k] += wt * v;
                        global_w[k] += wt;
                    }
                }
            }
        }
        ++region_nlearn[static_cast<std::size_t>(region[ui])];
        ++region_nlearn[static_cast<std::size_t>(ring[ui])];
        ++global_nlearn;
    }

    if (learned != nullptr) {
        learned->half = peak_half;
        learned->oversample = S;
        learned->n_bins_x = n_bins_x;
        learned->n_bins_y = n_bins_y;
        learned->region_sum = region_sum;
        learned->region_w = region_w;
        learned->region_n = region_nlearn;
        learned->global_sum = global_sum;
        learned->global_w = global_w;
        learned->global_n = global_nlearn;
    }
    // Pool learners from earlier frames (same peak size, sampling, tiling).
    if (const auto& prior = settings.prior;
        prior && !prior->empty() && prior->half == peak_half && prior->oversample == S &&
        prior->n_bins_x == n_bins_x && prior->n_bins_y == n_bins_y &&
        prior->region_sum.size() == region_sum.size()) {
        for (std::size_t k = 0; k < region_sum.size(); ++k) {
            region_sum[k] += prior->region_sum[k];
            region_w[k] += prior->region_w[k];
        }
        for (int r = 0; r < n_slots; ++r) {
            region_nlearn[static_cast<std::size_t>(r)] += prior->region_n[static_cast<std::size_t>(r)];
        }
        for (std::size_t k = 0; k < nf; ++k) {
            global_sum[k] += prior->global_sum[k];
            global_w[k] += prior->global_w[k];
        }
        global_nlearn += prior->global_n;
    }

    // Oversampling needs ~S² more learners than an integer-pixel profile to
    // fill every cell; regions short of that use the global profile, and a
    // Gaussian seed stands in until enough learners exist at all.
    const int min_fine = std::max(settings.min_learners, settings.min_learners * S * S / 4);
    std::vector<double> seed_sum;
    std::vector<double> seed_w;
    if (global_nlearn < min_fine) {
        seed_sum.assign(nf, 0.0);
        seed_w.assign(nf, 1.0);
        const double sg = std::max(settings.fallback_sigma_px, 0.3);
        for (int j = 0; j < grid.side; ++j) {
            for (int ci = 0; ci < grid.side; ++ci) {
                const double x = (ci + 0.5) / S - peak_half - 0.5;
                const double y = (j + 0.5) / S - peak_half - 0.5;
                seed_sum[static_cast<std::size_t>(j * grid.side + ci)] =
                    std::exp(-(x * x + y * y) / (2.0 * sg * sg)) / (2.0 * 3.14159265358979323846 * sg * sg);
            }
        }
    }
    std::vector<double> P(static_cast<std::size_t>(np));

    for (int i = 0; i < n; ++i) {
        const auto ui = static_cast<std::size_t>(i);
        ProfileObservation obs;
        static_cast<MeasObservation&>(obs) = meas[ui];
        obs.summation_intensity = obs.intensity;
        obs.summation_sigma = obs.sigma;
        obs.profile_region = region[ui];

        if (!obs.ok) {
            out.push_back(obs);
            continue;
        }
        const int c0 = col[ui];
        const int r0 = row[ui];
        const int reg = region[ui];
        const std::vector<double>* psum = &global_sum;
        const std::vector<double>* pw = &global_w;
        std::size_t off = 0;
        if (global_nlearn < min_fine) {
            psum = &seed_sum;
            pw = &seed_w;
        } else if (region_nlearn[static_cast<std::size_t>(reg)] >= min_fine) {
            psum = &region_sum;
            pw = &region_w;
            off = static_cast<std::size_t>(reg) * nf;
        } else if (region_nlearn[static_cast<std::size_t>(ring[ui])] >= min_fine) {
            psum = &region_sum;
            pw = &region_w;
            off = static_cast<std::size_t>(ring[ui]) * nf;
        }
        const double* rs = psum->data() + off;
        const double* rw = pw->data() + off;
        double pr3 = 0.0;
        double pr4 = 0.0;
        for (int dr = -peak_half; dr <= peak_half; ++dr) {
            for (int dc = -peak_half; dc <= peak_half; ++dc) {
                const double v = grid.sample(rs, rw, dc - fx[ui], dr - fy[ui]);
                P[static_cast<std::size_t>((dr + peak_half) * side + (dc + peak_half))] = v;
                pr3 += v;
                pr4 += v * v;
            }
        }

        const BackgroundPlane bg = background(i);
        obs.bg_a = bg.a;
        obs.bg_b = bg.b;
        obs.bg_c = bg.c;
        if (!bg.ok || pr4 < 1e-18) {
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
                const std::int32_t raw = pixel_at(image, c0 + dc, r0 + dr);
                const double b = bg.a * static_cast<double>(dc) +
                                 bg.b * static_cast<double>(dr) + bg.c;
                // Background plane under the peak (for Poisson σ model).
                bg_under += b;
                ++n_peak;
                if (raw < 0) {
                    continue;
                }
                const double Pk = P[static_cast<std::size_t>((dr + peak_half) * side + (dc + peak_half))];
                numer += Pk * (static_cast<double>(raw) - b);
                pr4_valid += Pk * Pk;
                p_valid += Pk;
                p2b_valid += Pk * Pk * std::max(b, 0.0);
                p3_valid += Pk * Pk * Pk;
            }
        }
        if (pr4_valid < 1e-18) {
            obs.profile_ok = false;
            out.push_back(obs);
            continue;
        }
        // I = k·ΣP: the volume inside the box of a profile normalised to the
        // whole spot (≈1 unless the box truncates it at this offset).
        const double iprf = numer * pr3 / pr4_valid;
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
        double var = pr3 * pr3 * (var_k_peak + var_k_bg);
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
