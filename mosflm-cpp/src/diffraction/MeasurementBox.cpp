#include "mosflm/diffraction/MeasurementBox.hpp"

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
    // Peak footprint ≈ half SEPARATION (mm → px), at least 2 px.
    const int hx = static_cast<int>(std::lround(0.5 * box.separation_x_mm / px));
    const int hy = static_cast<int>(std::lround(0.5 * box.separation_y_mm / py));
    return std::max(2, std::min(hx, hy));
}

}  // namespace (helpers above stay local)

void refine_centroid_com(
    const image::DiffractionImage& image,
    int& col,
    int& row,
    int half,
    double max_shift_px
) {
    if (half < 1 || max_shift_px <= 0.0) {
        return;
    }
    const int w = image.metadata.width;
    const int h = image.metadata.height;
    double sum = 0.0;
    double sx = 0.0;
    double sy = 0.0;
    double vmax = 0.0;
    int n_pix = 0;
    for (int dr = -half; dr <= half; ++dr) {
        for (int dc = -half; dc <= half; ++dc) {
            const int c = col + dc;
            const int r = row + dr;
            if (c < 0 || r < 0 || c >= w || r >= h) {
                continue;
            }
            const auto idx =
                static_cast<std::size_t>(r) * static_cast<std::size_t>(w) +
                static_cast<std::size_t>(c);
            const double v = static_cast<double>(image.pixels[idx]);
            if (v <= 0.0) {
                continue;
            }
            sum += v;
            sx += v * static_cast<double>(c);
            sy += v * static_cast<double>(r);
            vmax = std::max(vmax, v);
            ++n_pix;
        }
    }
    if (sum < 1e-6 || n_pix < 4) {
        return;
    }
    // Skip weak / flat windows — COM walks to noise and hurts I correlation.
    const double mean = sum / static_cast<double>(n_pix);
    if (!(vmax > 1.5 * mean + 1.0)) {
        return;
    }
    const double c_new = sx / sum;
    const double r_new = sy / sum;
    const double dx = c_new - static_cast<double>(col);
    const double dy = r_new - static_cast<double>(row);
    const double shift = std::hypot(dx, dy);
    if (shift > max_shift_px) {
        const double s = max_shift_px / shift;
        col = static_cast<int>(std::lround(static_cast<double>(col) + dx * s));
        row = static_cast<int>(std::lround(static_cast<double>(row) + dy * s));
    } else {
        col = static_cast<int>(std::lround(c_new));
        row = static_cast<int>(std::lround(r_new));
    }
}

namespace {

bool in_neighbour_peak(
    int c,
    int r,
    int peak_half,
    const std::vector<int>& neigh_cols,
    const std::vector<int>& neigh_rows
) {
    for (std::size_t k = 0; k < neigh_cols.size(); ++k) {
        if (std::abs(c - neigh_cols[k]) <= peak_half &&
            std::abs(r - neigh_rows[k]) <= peak_half) {
            return true;
        }
    }
    return false;
}

struct BgSample {
    double p{};
    double q{};
    double od{};
};

std::vector<BgSample> collect_bg_samples(
    const image::DiffractionImage& image,
    int c0,
    int r0,
    int peak_half,
    int outer,
    const std::vector<int>& neigh_cols,
    const std::vector<int>& neigh_rows,
    const std::vector<std::uint8_t>* trusted_mask = nullptr,
    int trusted_width = 0
) {
    std::vector<BgSample> samples;
    samples.reserve(static_cast<std::size_t>((2 * outer + 1) * (2 * outer + 1)));
    const int w = image.metadata.width;
    for (int dr = -outer; dr <= outer; ++dr) {
        for (int dc = -outer; dc <= outer; ++dc) {
            if (std::abs(dr) <= peak_half && std::abs(dc) <= peak_half) {
                continue;
            }
            const int c = c0 + dc;
            const int r = r0 + dr;
            if (c < 0 || r < 0 || c >= w || r >= image.metadata.height) {
                continue;
            }
            if (trusted_mask != nullptr && trusted_width > 0 &&
                !trusted_mask->empty()) {
                const auto idx =
                    static_cast<std::size_t>(r) *
                        static_cast<std::size_t>(trusted_width) +
                    static_cast<std::size_t>(c);
                if (idx >= trusted_mask->size() || (*trusted_mask)[idx] == 0) {
                    continue;
                }
            }
            if (in_neighbour_peak(c, r, peak_half, neigh_cols, neigh_rows)) {
                continue;
            }
            samples.push_back(
                {static_cast<double>(dc), static_cast<double>(dr),
                 static_cast<double>(pixel_at(image, c, r))}
            );
        }
    }
    return samples;
}

BackgroundPlane solve_plane(const std::vector<BgSample>& samples) {
    BackgroundPlane out;
    out.n_used = static_cast<int>(samples.size());
    if (samples.empty()) {
        return out;
    }
    double spp = 0, sqq = 0, spq = 0, sp = 0, sq = 0, sn = 0;
    double spod = 0, sqod = 0, sod = 0;
    for (const auto& s : samples) {
        spp += s.p * s.p;
        sqq += s.q * s.q;
        spq += s.p * s.q;
        sp += s.p;
        sq += s.q;
        sn += 1.0;
        spod += s.p * s.od;
        sqod += s.q * s.od;
        sod += s.od;
    }
    if (sn < 6.0) {
        out.c = sod / sn;
        out.ok = true;
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

}  // namespace

BackgroundPlane fit_background_plane(
    const image::DiffractionImage& image,
    int c0,
    int r0,
    int peak_half,
    int outer,
    const std::vector<int>& neigh_cols,
    const std::vector<int>& neigh_rows,
    double reject_nsigma,
    const std::vector<std::uint8_t>* trusted_mask,
    int trusted_width
) {
    auto samples = collect_bg_samples(
        image, c0, r0, peak_half, outer, neigh_cols, neigh_rows, trusted_mask,
        trusted_width
    );
    const int n0 = static_cast<int>(samples.size());
    BackgroundPlane plane = solve_plane(samples);
    if (!plane.ok) {
        plane.n_rejected = 0;
        return plane;
    }

    if (reject_nsigma > 0.0 && samples.size() >= 6) {
        for (int pass = 0; pass < 3; ++pass) {
            std::vector<double> resid;
            resid.reserve(samples.size());
            for (const auto& s : samples) {
                resid.push_back(s.od - plane.at(s.p, s.q));
            }
            double mean = 0.0;
            for (double v : resid) {
                mean += v;
            }
            mean /= static_cast<double>(resid.size());
            double var = 0.0;
            for (double v : resid) {
                const double d = v - mean;
                var += d * d;
            }
            const double sd =
                std::sqrt(var / static_cast<double>(resid.size()));
            if (sd < 1e-9) {
                break;
            }
            const double lim = reject_nsigma * sd;
            std::vector<BgSample> kept;
            kept.reserve(samples.size());
            for (std::size_t i = 0; i < samples.size(); ++i) {
                if (std::fabs(resid[i] - mean) <= lim) {
                    kept.push_back(samples[i]);
                }
            }
            if (kept.size() == samples.size() || kept.size() < 6) {
                samples.swap(kept);
                break;
            }
            samples.swap(kept);
            plane = solve_plane(samples);
            if (!plane.ok) {
                break;
            }
        }
    }

    plane.n_rejected = n0 - static_cast<int>(samples.size());
    plane.n_used = static_cast<int>(samples.size());
    return plane;
}

std::vector<MeasObservation> MeasIntegrator::integrate(
    const IntegrationRequest& request,
    const MeasurementBoxSettings& box
) const {
    const auto& image = request.image;
    const auto& detector = request.detector;
    const auto& predictions = request.reflections;
    const auto& settings = request.settings;

    std::vector<MeasObservation> out;
    const int n = static_cast<int>(predictions.size());
    if (n == 0) {
        return out;
    }

    const int peak_half = derive_peak_half(detector, box);
    const int rim = std::max(0, box.bg_rim_px);
    const int outer = peak_half + rim;
    if (peak_half < 0) {
        return out;
    }

    std::vector<int> col(static_cast<std::size_t>(n));
    std::vector<int> row(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        to_pixel(
            detector, settings, predictions[static_cast<std::size_t>(i)].x_mm,
            predictions[static_cast<std::size_t>(i)].y_mm,
            col[static_cast<std::size_t>(i)], row[static_cast<std::size_t>(i)]
        );
        if (box.refine_centroid) {
            refine_centroid_com(
                image, col[static_cast<std::size_t>(i)],
                row[static_cast<std::size_t>(i)], box.centroid_half_px,
                box.centroid_max_shift_px
            );
        }
    }

    const double sep_x_px =
        box.separation_x_mm / std::max(detector.pixel_size_x_mm, 1e-9);
    const double sep_y_px =
        box.separation_y_mm / std::max(detector.pixel_size_y_mm, 1e-9);
    std::vector<char> overlapped(static_cast<std::size_t>(n), 0);
    std::vector<std::vector<int>> neighbours(static_cast<std::size_t>(n));
    if (box.mark_overlaps) {
        const double neigh_x =
            std::max(sep_x_px, static_cast<double>(2 * outer + 1));
        const double neigh_y =
            std::max(sep_y_px, static_cast<double>(2 * outer + 1));
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                const double dx = static_cast<double>(
                    col[static_cast<std::size_t>(i)] -
                    col[static_cast<std::size_t>(j)]
                );
                const double dy = static_cast<double>(
                    row[static_cast<std::size_t>(i)] -
                    row[static_cast<std::size_t>(j)]
                );
                if (std::fabs(dx) < sep_x_px && std::fabs(dy) < sep_y_px) {
                    overlapped[static_cast<std::size_t>(i)] = 1;
                    overlapped[static_cast<std::size_t>(j)] = 1;
                }
                if (std::fabs(dx) < neigh_x && std::fabs(dy) < neigh_y) {
                    neighbours[static_cast<std::size_t>(i)].push_back(j);
                    neighbours[static_cast<std::size_t>(j)].push_back(i);
                }
            }
        }
    }

    out.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        MeasObservation r;
        const auto& pred = predictions[static_cast<std::size_t>(i)];
        r.h = pred.h;
        r.k = pred.k;
        r.l = pred.l;
        r.x_mm = pred.x_mm;
        r.y_mm = pred.y_mm;
        r.overlapped = overlapped[static_cast<std::size_t>(i)] != 0;
        r.box_half_px = outer;

        if (r.overlapped && !box.integrate_overlaps) {
            r.ok = false;
            out.push_back(r);
            continue;
        }

        const int c0 = col[static_cast<std::size_t>(i)];
        const int r0 = row[static_cast<std::size_t>(i)];
        std::vector<int> nc;
        std::vector<int> nr;
        for (int j : neighbours[static_cast<std::size_t>(i)]) {
            nc.push_back(col[static_cast<std::size_t>(j)]);
            nr.push_back(row[static_cast<std::size_t>(j)]);
        }

        const std::vector<std::uint8_t>* tmask =
            box.trusted_mask.empty() ? nullptr : &box.trusted_mask;
        const int tw = box.trusted_width;
        BackgroundPlane bg;
        if (box.use_bg_plane) {
            bg = fit_background_plane(
                image, c0, r0, peak_half, outer, nc, nr, box.bg_reject_nsigma,
                tmask, tw
            );
        } else {
            // Flat mean of rim (legacy path).
            auto samples = collect_bg_samples(
                image, c0, r0, peak_half, outer, nc, nr, tmask, tw
            );
            if (box.bg_reject_nsigma > 0.0 && samples.size() >= 4) {
                for (int pass = 0; pass < 3; ++pass) {
                    double sum = 0.0;
                    for (const auto& s : samples) {
                        sum += s.od;
                    }
                    const double mean = sum / static_cast<double>(samples.size());
                    double var = 0.0;
                    for (const auto& s : samples) {
                        const double d = s.od - mean;
                        var += d * d;
                    }
                    const double sd =
                        std::sqrt(var / static_cast<double>(samples.size()));
                    if (sd < 1e-9) {
                        break;
                    }
                    const double lim = box.bg_reject_nsigma * sd;
                    std::vector<BgSample> kept;
                    for (const auto& s : samples) {
                        if (std::fabs(s.od - mean) <= lim) {
                            kept.push_back(s);
                        }
                    }
                    if (kept.size() == samples.size() || kept.size() < 4) {
                        samples.swap(kept);
                        break;
                    }
                    samples.swap(kept);
                }
            }
            bg.n_used = static_cast<int>(samples.size());
            if (!samples.empty()) {
                double sum = 0.0;
                for (const auto& s : samples) {
                    sum += s.od;
                }
                bg.c = sum / static_cast<double>(samples.size());
                bg.ok = true;
            }
        }

        r.n_bg_rejected = bg.n_rejected;
        r.npix_bg = bg.n_used;
        r.bg_a = bg.a;
        r.bg_b = bg.b;
        r.bg_c = bg.c;

        double peak_sum = 0.0;
        double bg_under_peak = 0.0;
        int n_peak = 0;
        for (int dr = -peak_half; dr <= peak_half; ++dr) {
            for (int dc = -peak_half; dc <= peak_half; ++dc) {
                const int c = c0 + dc;
                const int rr = r0 + dr;
                if (c < 0 || rr < 0 || c >= image.metadata.width ||
                    rr >= image.metadata.height) {
                    continue;
                }
                peak_sum += static_cast<double>(pixel_at(image, c, rr));
                if (bg.ok) {
                    bg_under_peak += bg.at(
                        static_cast<double>(dc), static_cast<double>(dr)
                    );
                }
                ++n_peak;
            }
        }
        r.npix_peak = n_peak;
        if (n_peak == 0) {
            r.ok = false;
            out.push_back(r);
            continue;
        }

        if (!bg.ok) {
            bg_under_peak = 0.0;
        }
        r.intensity = peak_sum - bg_under_peak;
        // Counting-stats stand-in: peak variance + bg contribution.
        const double var =
            std::fabs(peak_sum) +
            (r.npix_bg > 0
                 ? (static_cast<double>(n_peak * n_peak) /
                    static_cast<double>(r.npix_bg)) *
                       std::fabs(bg_under_peak / std::max(n_peak, 1))
                 : 0.0);
        r.sigma = std::sqrt(std::max(var, 1.0));
        if (settings.apply_lp) {
            apply_lp_inplace(
                detector, r.x_mm, r.y_mm, r.intensity, r.sigma,
                settings.polarization_fraction
            );
        }
        r.ok = true;
        out.push_back(r);
    }
    return out;
}

}  // namespace mosflm::diffraction
