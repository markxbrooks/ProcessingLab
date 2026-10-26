/**
 * 3D shoebox extraction and Rossmann profile fit.
 */

#include "mosflm/diffraction/Profile3d.hpp"
#include "mosflm/diffraction/Shoebox.hpp"
#include "test_util.hpp"

#include <cmath>
#include <cstdint>
#include <random>
#include <string>
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
    // σ_M / osc is the φ profile width. This σ_M keeps that width at
    // 1/2.355 frames, where the empty-box χ² below was calibrated.
    s3.mosaic_spread_deg = 0.15 / 2.355;
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

    // Sparse Eiger background (λ = 0.1 counts/voxel, no peak): I must be
    // unbiased with ~half negative, and σ must match the scatter (χ² ≈ 1).
    // A symmetric bg clip zeroes the bg here and gives I ≈ +24, never < 0.
    {
        Profile3dSettings sp = s3;
        sp.use_phi_bg = false;
        sp.box.bg_reject_nsigma = 3.0;
        sp.min_learners = 1000000;  // fixed Gaussian profile
        std::mt19937 rng(12345);
        std::poisson_distribution<int> pois(0.1);
        std::vector<Shoebox> boxes;
        for (int n = 0; n < 2000; ++n) {
            auto box = make_empty_box(3, 2, 5);
            for (auto& v : box.pixels) {
                v = pois(rng);
            }
            boxes.push_back(std::move(box));
        }
        const auto obs = Profile3dFitter{}.fit(boxes, det, iset, sp);
        double sum_i = 0.0, chi2_sum = 0.0, chi2_prf = 0.0;
        int n_ok = 0, n_neg = 0, n_zero = 0;
        for (const auto& o : obs) {
            if (!o.ok || !o.profile_ok) {
                continue;
            }
            ++n_ok;
            sum_i += o.summation_intensity;
            n_neg += o.summation_intensity < 0.0 ? 1 : 0;
            n_zero += o.summation_intensity == 0.0 ? 1 : 0;
            chi2_sum += std::pow(o.summation_intensity / o.summation_sigma, 2);
            chi2_prf += std::pow(o.profile_intensity / o.profile_sigma, 2);
        }
        CHECK(n_ok > 1900);
        const double mean_i = sum_i / n_ok;
        CHECK_MSG(std::fabs(mean_i) < 0.6, "sparse bg mean I " + mosflm_test::fmt(mean_i));
        CHECK_MSG(
            n_neg > n_ok / 3, "sparse bg negatives " + std::to_string(n_neg)
        );
        CHECK_MSG(n_zero < n_ok / 20, "sparse bg zeros " + std::to_string(n_zero));
        chi2_sum /= n_ok;
        chi2_prf /= n_ok;
        CHECK_MSG(
            chi2_sum > 0.8 && chi2_sum < 1.25,
            "summation chi2 " + mosflm_test::fmt(chi2_sum)
        );
        CHECK_MSG(
            chi2_prf > 0.8 && chi2_prf < 1.25,
            "profile chi2 " + mosflm_test::fmt(chi2_prf)
        );
    }

    // Strong spot on a moderate background: σ matches the repeat scatter.
    {
        Profile3dSettings sp = s3;
        sp.use_phi_bg = false;
        sp.box.bg_reject_nsigma = 3.0;
        sp.min_learners = 1000000;
        const int peak_half = 3;
        const int nphi = 5;
        const int cphi = nphi / 2;
        double psum = 0.0;
        for (int iphi = 0; iphi < nphi; ++iphi) {
            for (int dr = -peak_half; dr <= peak_half; ++dr) {
                for (int dc = -peak_half; dc <= peak_half; ++dc) {
                    psum += gaussian3(dc, dr, iphi - cphi, 1.2, 0.8);
                }
            }
        }
        std::mt19937 rng(777);
        std::vector<Shoebox> boxes;
        for (int n = 0; n < 1000; ++n) {
            auto box = make_empty_box(peak_half, 2, nphi);
            const int outer = box.outer_half;
            for (int iphi = 0; iphi < nphi; ++iphi) {
                for (int dr = -outer; dr <= outer; ++dr) {
                    for (int dc = -outer; dc <= outer; ++dc) {
                        double m = 2.0;
                        if (std::abs(dc) <= peak_half && std::abs(dr) <= peak_half) {
                            m += 2000.0 *
                                 gaussian3(dc, dr, iphi - cphi, 1.2, 0.8) / psum;
                        }
                        box.pixels[box.index(iphi, dr, dc)] =
                            std::poisson_distribution<int>(m)(rng);
                    }
                }
            }
            boxes.push_back(std::move(box));
        }
        const auto obs = Profile3dFitter{}.fit(boxes, det, iset, sp);
        std::vector<double> is, ip, ss, sprf;
        for (const auto& o : obs) {
            if (o.ok && o.profile_ok) {
                is.push_back(o.summation_intensity);
                ss.push_back(o.summation_sigma);
                ip.push_back(o.profile_intensity);
                sprf.push_back(o.profile_sigma);
            }
        }
        auto ratio = [](const std::vector<double>& v, const std::vector<double>& s) {
            double mean = 0.0, m2 = 0.0, ms = 0.0;
            for (double x : v) {
                mean += x;
            }
            mean /= static_cast<double>(v.size());
            for (std::size_t i = 0; i < v.size(); ++i) {
                m2 += (v[i] - mean) * (v[i] - mean);
                ms += s[i];
            }
            const double sd = std::sqrt(m2 / static_cast<double>(v.size() - 1));
            return (ms / static_cast<double>(s.size())) / sd;
        };
        CHECK(is.size() > 950);
        const double rs = ratio(is, ss);
        const double rp = ratio(ip, sprf);
        CHECK_MSG(rs > 0.85 && rs < 1.15, "strong summation σ/sd " + mosflm_test::fmt(rs));
        CHECK_MSG(rp > 0.85 && rp < 1.15, "strong profile σ/sd " + mosflm_test::fmt(rp));
    }

    return mosflm_test::finish("profile3d");
}
