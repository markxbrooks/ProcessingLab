/**
 * sum_partial_slices: rocking-curve slices → one reflection.
 */
#include "mosflm/diffraction/PartialSummation.hpp"
#include "test_util.hpp"

#include <cmath>

int main() {
    using mosflm::crystal::MosflmReflection;
    using mosflm::diffraction::sum_partial_slices;
    using mosflm_test::close;

    auto slice = [](int batch, float frac, float I, float phi, int flag) {
        MosflmReflection r;
        r.h = r.h0 = 1;
        r.k = r.k0 = 2;
        r.l = r.l0 = 3;
        r.batch = batch;
        r.rot = 10.0f;
        r.fractioncalc = frac;
        r.intensity = r.ipr = I;
        r.sigma = r.sigipr = 2.0f;
        r.phi_obs = phi;
        r.flag = flag;
        r.xdet_obs = 100.0f + static_cast<float>(batch);
        r.ydet_obs = 200.0f;
        return r;
    };
    std::vector<MosflmReflection> rows{
        slice(66, 0.25f, 25.0f, 9.925f, 0), slice(67, 0.5f, 50.0f, 10.075f, 0),
        slice(68, 0.24f, 24.0f, 10.225f, 4),  // weak-tail neg-profile bit dropped
        slice(68, 0.05f, -1.0f, 10.225f, 16),  // tiny tail: overlap bit ignored
    };
    // Same HKL 180° away is a different reflection.
    auto other = slice(1266, 1.0f, 7.0f, 190.0f, 0);
    other.rot = 190.0f;
    rows.push_back(other);

    const auto out = sum_partial_slices(rows);
    CHECK(out.size() == 2u);
    const auto& r = out[0];
    CHECK(r.batch == 67);
    CHECK(close(r.intensity, 98.0, 1e-4));
    CHECK(close(r.sigma, 4.0, 1e-4));
    CHECK(close(r.fractioncalc, 1.0, 1e-6));  // 1.04 capped
    CHECK(r.mpart == 0);
    CHECK(r.flag == 0);

    // Bit 4 is recomputed on the sum: within 5σ it stays clear, below −5σ
    // it is set. Each slice σ is 2, so two slices give σ = √8.
    {
        auto mild_a = slice(1, 0.5f, -2.0f, 10.0f, 4);
        auto mild_b = slice(2, 0.5f, -2.0f, 10.2f, 0);
        const auto mild = sum_partial_slices({mild_a, mild_b});
        CHECK(mild.size() == 1u);
        CHECK(mild[0].flag == 0);
        auto hard_a = slice(1, 0.5f, -40.0f, 10.0f, 0);
        auto hard_b = slice(2, 0.5f, -40.0f, 10.2f, 0);
        const auto hard = sum_partial_slices({hard_a, hard_b});
        CHECK(hard.size() == 1u);
        CHECK((hard[0].flag & 4) != 0);
    }
    // Re-centred with signed weights in the window: the −1 tail slice counts.
    const double phi = (25 * 9.925 + 50 * 10.075 + 23 * 10.225) / 98.0;
    CHECK(close(r.phi_obs, phi, 1e-4));
    CHECK(close(r.xdet_obs, (25 * 166 + 50 * 167 + 24 * 168) / 99.0, 1e-3));
    CHECK(out[1].batch == 1266 && close(out[1].intensity, 7.0, 1e-6));

    // φ spread uses signed weights: ± noise slices far out in the reach
    // cancel instead of inflating it (clipped weights fed a MOSAIC runaway).
    std::vector<MosflmReflection> rock{
        slice(66, 0.25f, 25.0f, 9.925f, 0), slice(67, 0.5f, 50.0f, 10.075f, 0),
        slice(68, 0.25f, 25.0f, 10.225f, 0),
        slice(56, 0.0f, 6.0f, 8.425f, 0), slice(57, 0.0f, -6.0f, 8.575f, 0),
        slice(78, 0.0f, 6.0f, 11.725f, 0), slice(77, 0.0f, -6.0f, 11.575f, 0),
    };
    const auto rs = sum_partial_slices(rock);
    CHECK(rs.size() == 1u);
    // Signed: (1.125 + 2·6·(1.65² − 1.5²)) / 100 → 0.261°; clipped
    // weights would give √((1.125 + 12·1.65²) / 112) = 0.550°.
    CHECK(close(rs[0].phi_obs_sd, 0.2608, 1e-3));

    // Rocking-width fit: slices predicted at W = 0.5° (reach 4σ + 2 frames)
    // of curves that are really W = 1.0° (ζ = 0.5) — the reach truncates
    // them at ~2.5σ, which the fit models instead of reading as narrow.
    {
        using mosflm::diffraction::fit_rocking_width;
        const double osc = 0.2, zeta = 0.5, w_now = 0.5, w_true = 1.0;
        const double s_now = w_now / (2.355 * zeta);
        const double s_true = w_true / (2.355 * zeta);
        auto cdf = [](double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); };
        std::vector<MosflmReflection> rows2;
        for (int n = 0; n < 60; ++n) {
            const double phic = 10.0 + 2.0 * n + 0.037 * n;  // off-grid φcalc
            const double reach = 4.0 * s_now + 2.0 * osc;
            for (int f = static_cast<int>(std::floor((phic - reach) / osc));
                 f * osc <= phic + reach; ++f) {
                const double f0 = f * osc - phic, f1 = f0 + osc;
                auto r = slice(f, 0.0f, 0.0f, 0.0f, 0);
                r.h = r.h0 = n;
                r.rot = static_cast<float>(phic);
                r.width = static_cast<float>(6.0 * s_now);
                r.fractioncalc = static_cast<float>(cdf(f1 / s_now) - cdf(f0 / s_now));
                r.intensity = static_cast<float>(
                    2000.0 * (cdf(f1 / s_true) - cdf(f0 / s_true)) + ((f & 1) ? 1.5 : -1.5)
                );
                r.phi_obs = static_cast<float>(f * osc + 0.5 * osc);
                rows2.push_back(r);
            }
        }
        const auto fit = fit_rocking_width(rows2, w_now, osc);
        CHECK(fit.has_value() && fit->n_reflections == 60u);
        CHECK(fit.has_value() && close(fit->width_deg, w_true, 0.03));
        // Too few strong reflections → no fit.
        CHECK(!fit_rocking_width({rows2.begin(), rows2.begin() + 30}, w_now, osc).has_value());
    }

    // Predicted too wide (W = 2.0°, so many tail frames) with a small
    // positive background bias on every slice: the bias must not read as
    // heavier tails. Without a baseline in the fit this fed back through
    // MOSAIC → reach and ran away (76c1b: 0.96° → 2.0°; XDS 1.04°).
    {
        using mosflm::diffraction::fit_rocking_width;
        const double osc = 0.2, zeta = 0.5, w_now = 2.0, w_true = 1.0;
        const double s_now = w_now / (2.355 * zeta);
        const double s_true = w_true / (2.355 * zeta);
        auto cdf = [](double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); };
        std::vector<MosflmReflection> rows3;
        for (int n = 0; n < 60; ++n) {
            const double phic = 10.0 + 4.0 * n + 0.037 * n;
            const double reach = 4.0 * s_now + 2.0 * osc;
            for (int f = static_cast<int>(std::floor((phic - reach) / osc)); f * osc <= phic + reach; ++f) {
                const double f0 = f * osc - phic, f1 = f0 + osc;
                auto r = slice(f, 0.0f, 0.0f, 0.0f, 0);
                r.h = r.h0 = n;
                r.rot = static_cast<float>(phic);
                r.width = static_cast<float>(6.0 * s_now);
                r.fractioncalc = static_cast<float>(cdf(f1 / s_now) - cdf(f0 / s_now));
                r.intensity = static_cast<float>(2000.0 * (cdf(f1 / s_true) - cdf(f0 / s_true)) + 8.0);
                r.phi_obs = static_cast<float>(f * osc + 0.5 * osc);
                rows3.push_back(r);
            }
        }
        const auto fit = fit_rocking_width(rows3, w_now, osc);
        CHECK(fit.has_value() && close(fit->width_deg, w_true, 0.05));
    }

    // Heavy-tailed rocking curves (a Gaussian core with 20% of the flux in a
    // broad Lorentzian, as on 76c1b): the width must not depend on how far
    // the prediction reached. A Gaussian fitted to every slice grew with the
    // width predicted with, so MOSAIC crept up at every re-estimate.
    {
        using mosflm::diffraction::fit_rocking_width;
        const double osc = 0.2, zeta = 0.6, w_core = 0.8, hwhm_tail = 1.5;
        auto cdf = [](double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); };
        const auto lorentz_cdf = [&](double x) { return 0.5 + std::atan(x * zeta / hwhm_tail) / 3.14159265358979; };
        const auto fitted_with = [&](double w_now) {
            const double s_now = w_now / (2.355 * zeta);
            const double s_core = w_core / (2.355 * zeta);
            std::vector<MosflmReflection> rows4;
            for (int n = 0; n < 80; ++n) {
                const double phic = 10.0 + 7.0 * n + 0.037 * n;
                const double reach = 4.0 * s_now + 2.0 * osc;
                for (int f = static_cast<int>(std::floor((phic - reach) / osc)); f * osc <= phic + reach; ++f) {
                    const double f0 = f * osc - phic, f1 = f0 + osc;
                    auto r = slice(f, 0.0f, 0.0f, 0.0f, 0);
                    r.h = r.h0 = n;
                    r.rot = static_cast<float>(phic);
                    r.width = static_cast<float>(6.0 * s_now);
                    r.fractioncalc = static_cast<float>(cdf(f1 / s_now) - cdf(f0 / s_now));
                    r.intensity = static_cast<float>(
                        4000.0 * (0.8 * (cdf(f1 / s_core) - cdf(f0 / s_core)) +
                                  0.2 * (lorentz_cdf(f1) - lorentz_cdf(f0))));
                    r.phi_obs = static_cast<float>(f * osc + 0.5 * osc);
                    rows4.push_back(r);
                }
            }
            const auto fit = fit_rocking_width(rows4, w_now, osc);
            return fit ? fit->width_deg : -1.0;
        };
        const double narrow = fitted_with(0.6);
        const double wide = fitted_with(1.6);
        CHECK_MSG(narrow > 0.0 && wide > 0.0 && std::fabs(narrow - wide) < 0.1,
                  "heavy tails: " + std::to_string(narrow) + " vs " + std::to_string(wide));
        CHECK_MSG(std::fabs(wide - w_core) < 0.2, "core width " + std::to_string(wide));
    }

    // Incomplete edge reflection stays partial.
    const auto edge = sum_partial_slices({slice(1, 0.3f, 30.0f, 0.075f, 0)});
    CHECK(edge.size() == 1u && edge[0].mpart == 1 && close(edge[0].fractioncalc, 0.3, 1e-6));
    // Only a tail slice survived (core frames masked) → reflection dropped.
    const auto tail_only = sum_partial_slices({slice(70, 0.0005f, 0.0f, 10.5f, 0)});
    CHECK(tail_only.empty());
    return mosflm_test::finish("partial_summation");
}
