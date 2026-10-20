/**
 * 3D shoebox extraction and Rossmann profile fit.
 */

#include "mosflm/diffraction/Profile3d.hpp"
#include "mosflm/diffraction/Shoebox.hpp"
#include "test_util.hpp"

#include <cmath>
#include <cstdint>
#include <vector>

using mosflm::DetectorGeometry;
using mosflm::diffraction::FrameCache;
using mosflm::diffraction::IntegrationSettings;
using mosflm::diffraction::Profile3dFitter;
using mosflm::diffraction::Profile3dSettings;
using mosflm::diffraction::Shoebox;
using mosflm::diffraction::odd_nadd;
using mosflm::image::DiffractionImage;

namespace {

Shoebox make_empty_box(int peak_half, int bg_rim, int nphi) {
    Shoebox box;
    box.h = 1;
    box.k = 0;
    box.l = 0;
    box.col0 = 32;
    box.row0 = 32;
    box.peak_half = peak_half;
    box.outer_half = peak_half + bg_rim;
    box.nphi = nphi;
    box.peak_frame = 3;
    box.first_frame = 3 - nphi / 2;
    box.width = 64;
    box.height = 64;
    box.x_mm = 0.0;
    box.y_mm = 0.0;
    box.phi_calc_deg = 0.3;
    const int side = box.side();
    const std::size_t n =
        static_cast<std::size_t>(nphi) * static_cast<std::size_t>(side) *
        static_cast<std::size_t>(side);
    box.pixels.assign(n, 0);
    box.valid.assign(n, 1);
    box.ok = true;
    return box;
}

double gaussian3(
    int dc,
    int dr,
    int dphi,
    double sig_xy,
    double sig_phi
) {
    const double r2 = static_cast<double>(dc * dc + dr * dr);
    const double p2 = static_cast<double>(dphi * dphi);
    return std::exp(-0.5 * r2 / (sig_xy * sig_xy)) *
           std::exp(-0.5 * p2 / (sig_phi * sig_phi));
}

}  // namespace

int main() {
    CHECK_EQ(odd_nadd(5), 5);
    CHECK_EQ(odd_nadd(4), 5);
    CHECK_EQ(odd_nadd(1), 1);

    {
        FrameCache cache(2);
        DiffractionImage a;
        a.metadata.width = 4;
        a.metadata.height = 4;
        a.pixels.assign(16, 1);
        DiffractionImage b = a;
        b.pixels.assign(16, 2);
        DiffractionImage c = a;
        c.pixels.assign(16, 3);
        cache.put("mem", 1, a);
        cache.put("mem", 2, b);
        CHECK_EQ(static_cast<int>(cache.size()), 2);
        cache.put("mem", 3, c);
        CHECK_EQ(static_cast<int>(cache.size()), 2);
        auto got = cache.get("mem", 3);
        CHECK(got != nullptr);
        CHECK_EQ(got->pixels[0], 3);
    }

    DetectorGeometry det;
    det.nfast = 64;
    det.nslow = 64;
    det.pixel_size_x_mm = 0.1;
    det.pixel_size_y_mm = 0.1;
    det.distance_mm = 100.0;
    det.wavelength_A = 1.0;
    IntegrationSettings iset;
    Profile3dSettings s3;
    s3.nadd = 5;
    s3.min_learners = 1;
    s3.learn_min_isigi = 0.5;
    s3.n_bins_x = 1;
    s3.n_bins_y = 1;
    s3.fallback_sigma_xy_px = 1.2;
    s3.fallback_sigma_phi_frames = 0.8;
    s3.use_phi_bg = true;
    s3.box.peak_half_px = 3;
    s3.box.bg_rim_px = 2;
    s3.box.bg_reject_nsigma = 0.0;
    s3.mosaic_spread_deg = 0.15;
    s3.osc_deg = 0.15;

    // Synthetic 3D Gaussian + planar+φ background recovers I.
    {
        const int peak_half = 3;
        const int nphi = 5;
        const int cphi = nphi / 2;
        auto box = make_empty_box(peak_half, 2, nphi);
        double psum = 0.0;
        for (int iphi = 0; iphi < nphi; ++iphi) {
            for (int dr = -peak_half; dr <= peak_half; ++dr) {
                for (int dc = -peak_half; dc <= peak_half; ++dc) {
                    psum += gaussian3(dc, dr, iphi - cphi, 1.2, 0.8);
                }
            }
        }
        const double i_true = 8000.0;
        const int outer = box.outer_half;
        for (int iphi = 0; iphi < nphi; ++iphi) {
            const double dphi = static_cast<double>(iphi - cphi);
            for (int dr = -outer; dr <= outer; ++dr) {
                for (int dc = -outer; dc <= outer; ++dc) {
                    const double bg = 20.0 + 0.4 * static_cast<double>(dc) +
                                      0.2 * static_cast<double>(dr) + 3.0 * dphi;
                    double peak = 0.0;
                    if (std::abs(dc) <= peak_half && std::abs(dr) <= peak_half) {
                        peak = i_true *
                               gaussian3(dc, dr, iphi - cphi, 1.2, 0.8) / psum;
                    }
                    box.pixels[box.index(iphi, dr, dc)] =
                        static_cast<std::int32_t>(std::lround(bg + peak));
                }
            }
        }
        const auto obs = Profile3dFitter{}.fit({box}, det, iset, s3);
        CHECK_EQ(static_cast<int>(obs.size()), 1);
        CHECK(obs[0].ok);
        CHECK_MSG(
            std::fabs(obs[0].intensity - i_true) / i_true < 0.12,
            "3D Gaussian I recover " + mosflm_test::fmt(obs[0].intensity)
        );
    }

    // Two-frame 2D+sum ≠ 3D when background varies with φ and tails are
    // biased high in the peak box (2D includes negative tail I; 3D P≈0 there).
    {
        const int peak_half = 3;
        const int nphi = 5;
        const int cphi = nphi / 2;
        auto box = make_empty_box(peak_half, 2, nphi);
        const int outer = box.outer_half;
        const int side_peak = 2 * peak_half + 1;
        const int n_peak_xy = side_peak * side_peak;
        double psum = 0.0;
        for (int dr = -peak_half; dr <= peak_half; ++dr) {
            for (int dc = -peak_half; dc <= peak_half; ++dc) {
                psum += gaussian3(dc, dr, 0, 1.2, 0.8);
            }
        }
        const double i_true = 6000.0;
        for (int iphi = 0; iphi < nphi; ++iphi) {
            for (int dr = -outer; dr <= outer; ++dr) {
                for (int dc = -outer; dc <= outer; ++dc) {
                    const bool peak =
                        std::abs(dc) <= peak_half && std::abs(dr) <= peak_half;
                    double v = 0.0;
                    if (iphi == cphi) {
                        const double bg = 20.0;
                        double pk = 0.0;
                        if (peak) {
                            pk = i_true * gaussian3(dc, dr, 0, 1.2, 0.8) / psum;
                        }
                        v = bg + pk;
                    } else if (peak) {
                        v = 40.0;  // peak box below rim
                    } else {
                        v = 90.0;  // high rim → 2D bg overestimate
                    }
                    box.pixels[box.index(iphi, dr, dc)] =
                        static_cast<std::int32_t>(std::lround(v));
                }
            }
        }

        double i2d = 0.0;
        for (int iphi = 0; iphi < nphi; ++iphi) {
            double rim = 0.0;
            int nr = 0;
            for (int dr = -outer; dr <= outer; ++dr) {
                for (int dc = -outer; dc <= outer; ++dc) {
                    if (std::abs(dc) <= peak_half && std::abs(dr) <= peak_half) {
                        continue;
                    }
                    rim += static_cast<double>(box.at(iphi, dr, dc));
                    ++nr;
                }
            }
            const double bg = rim / static_cast<double>(std::max(1, nr));
            double islice = 0.0;
            for (int dr = -peak_half; dr <= peak_half; ++dr) {
                for (int dc = -peak_half; dc <= peak_half; ++dc) {
                    islice += static_cast<double>(box.at(iphi, dr, dc)) - bg;
                }
            }
            i2d += islice;
        }
        (void)n_peak_xy;
        const auto obs = Profile3dFitter{}.fit({box}, det, iset, s3);
        CHECK(obs[0].ok);
        CHECK_MSG(
            std::fabs(obs[0].intensity - i2d) > 200.0,
            "2D+sum should differ from 3D; 2D=" + mosflm_test::fmt(i2d) +
                " 3D=" + mosflm_test::fmt(obs[0].intensity)
        );
        CHECK_MSG(
            std::fabs(obs[0].intensity - i_true) < std::fabs(i2d - i_true),
            "3D closer to true I than 2D+sum"
        );
    }

    return mosflm_test::finish("profile3d");
}
