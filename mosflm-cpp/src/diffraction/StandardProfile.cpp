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
        mosflm::cpp_mm_to_mosflm_xdet(
            x_mm, y_mm, settings.image_beam_x_mm, settings.image_beam_y_mm,
            detector.pixel_size_x_mm, detector.pixel_size_y_mm, xdet, ydet
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

struct BgPlane {
    double a{0.0};
    double b{0.0};
    double c{0.0};
    bool ok{false};
};

BgPlane fit_bg_plane(
    const image::DiffractionImage& image,
    int c0,
    int r0,
    int peak_half,
    int outer,
    const std::vector<int>& neigh_cols,
    const std::vector<int>& neigh_rows
) {
    double spp = 0, sqq = 0, spq = 0, sp = 0, sq = 0, sn = 0;
    double spod = 0, sqod = 0, sod = 0;
    for (int dr = -outer; dr <= outer; ++dr) {
        for (int dc = -outer; dc <= outer; ++dc) {
            if (std::abs(dr) <= peak_half && std::abs(dc) <= peak_half) {
                continue;
            }
            const int c = c0 + dc;
            const int r = r0 + dr;
            if (c < 0 || r < 0 || c >= image.metadata.width ||
                r >= image.metadata.height) {
                continue;
            }
            bool neigh = false;
            for (std::size_t k = 0; k < neigh_cols.size(); ++k) {
                if (std::abs(c - neigh_cols[k]) <= peak_half &&
                    std::abs(r - neigh_rows[k]) <= peak_half) {
                    neigh = true;
                    break;
                }
            }
            if (neigh) {
                continue;
            }
            const double p = static_cast<double>(dc);
            const double q = static_cast<double>(dr);
            const double od = static_cast<double>(pixel_at(image, c, r));
            spp += p * p;
            sqq += q * q;
            spq += p * q;
            sp += p;
            sq += q;
            sn += 1.0;
            spod += p * od;
            sqod += q * od;
            sod += od;
        }
    }
    BgPlane out;
    if (sn < 6.0) {
        if (sn > 0.0) {
            out.c = sod / sn;
            out.ok = true;
        }
        return out;
    }
    const double det =
        spp * (sqq * sn - sq * sq) - spq * (spq * sn - sq * sp) +
        sp * (spq * sq - sqq * sp);
    if (std::fabs(det) < 1e-12) {
        out.c = sod / sn;
        out.ok = true;
        return out;
    }
    out.a = (spod * (sqq * sn - sq * sq) - spq * (sqod * sn - sq * sod) +
             sp * (sqod * sq - sqq * sod)) /
            det;
    out.b = (spp * (sqod * sn - sq * sod) - spod * (spq * sn - sq * sp) +
             sp * (spq * sod - sqod * sp)) /
            det;
    out.c = (spp * (sqq * sod - sqod * sq) - spq * (spq * sod - sqod * sp) +
             spod * (spq * sq - sqq * sp)) /
            det;
    out.ok = true;
    return out;
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
    double sum = 0.0;
    for (double v : profile) {
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
        if (!m.ok || m.sigma <= 0.0) {
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
        const BgPlane bg =
            fit_bg_plane(image, c0, r0, peak_half, outer, nc, nr);
        if (!bg.ok) {
            continue;
        }
        double peak_max = 0.0;
        std::vector<double> local(static_cast<std::size_t>(np), 0.0);
        for (int dr = -peak_half; dr <= peak_half; ++dr) {
            for (int dc = -peak_half; dc <= peak_half; ++dc) {
                const int c = c0 + dc;
                const int r = r0 + dr;
                const double od = static_cast<double>(pixel_at(image, c, r));
                const double b =
                    bg.a * static_cast<double>(dc) +
                    bg.b * static_cast<double>(dr) + bg.c;
                const double v = std::max(0.0, od - b);
                const std::size_t idx = static_cast<std::size_t>(
                    (dr + peak_half) * side + (dc + peak_half)
                );
                local[idx] = v;
                peak_max = std::max(peak_max, v);
            }
        }
        if (peak_max < 1.0) {
            continue;
        }
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
        const BgPlane bg =
            fit_bg_plane(image, c0, r0, peak_half, outer, nc, nr);
        obs.bg_a = bg.a;
        obs.bg_b = bg.b;
        obs.bg_c = bg.c;

        if (!bg.ok || pr.pr4 < 1e-18) {
            obs.profile_ok = false;
            out.push_back(obs);
            continue;
        }

        double sumprod = 0.0;
        double spkod = 0.0;
        double spkodsq = 0.0;
        for (int dr = -peak_half; dr <= peak_half; ++dr) {
            for (int dc = -peak_half; dc <= peak_half; ++dc) {
                const int c = c0 + dc;
                const int rr = r0 + dr;
                const double od = static_cast<double>(pixel_at(image, c, rr));
                const double P = profile[static_cast<std::size_t>(
                    (dr + peak_half) * side + (dc + peak_half)
                )];
                sumprod += P * od;
                spkod += od;
                spkodsq += od * od;
            }
        }

        const double numer =
            sumprod - pr.pr1 * bg.a - pr.pr2 * bg.b - pr.pr3 * bg.c;
        const double iprf = numer * pr.pr3 / pr.pr4;
        obs.profile_intensity = iprf;
        const double xj = (pr.pr3 > 0.0) ? (iprf / pr.pr3) : 0.0;
        double var = xj * xj * pr.pr4 + std::fabs(spkodsq);
        var = std::max(var, 1.0);
        obs.profile_sigma = std::sqrt(var);
        obs.pkratio =
            std::fabs(obs.summation_intensity - iprf) /
            std::max(std::fabs(iprf), 1.0);
        obs.profile_ok = true;

        if (settings.prefer_profile) {
            obs.intensity = obs.profile_intensity;
            obs.sigma = obs.profile_sigma;
        }
        if (iset.apply_lp) {
            apply_lp_inplace(
                detector, obs.x_mm, obs.y_mm, obs.intensity, obs.sigma
            );
            apply_lp_inplace(
                detector, obs.x_mm, obs.y_mm, obs.profile_intensity,
                obs.profile_sigma
            );
            apply_lp_inplace(
                detector, obs.x_mm, obs.y_mm, obs.summation_intensity,
                obs.summation_sigma
            );
        }
        out.push_back(obs);
    }
    return out;
}

}  // namespace mosflm::diffraction
