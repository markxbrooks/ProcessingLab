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
    // Bit 512 is a spatial overlap kept for the mosaicity sample. The sum
    // drops that slice, so its intensity and flag do not reach the MTZ.
    {
        auto keep = slice(1, 0.6f, 40.0f, 10.0f, 0);
        auto drop = slice(2, 0.4f, 80.0f, 10.2f, 512);
        const auto summed = sum_partial_slices({keep, drop});
        CHECK(summed.size() == 1u);
        CHECK(close(summed[0].intensity, 40.0, 1e-4));
        CHECK(summed[0].flag == 0);
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

    // σ_M from the core of each rocking curve, inside a fixed |ε| ≤ 1° cap.
    // ζ = 0.5, true σ_M = 0.40° so σφ = 0.80°. The same curves predicted at
    // two different σ_M (different slice widths) must return the same value.
    {
        using mosflm::diffraction::fit_rocking_width;
        const double osc = 0.2;
        const double zeta = 0.5;
        const double sigma_true = 0.40;
        const double s_true = sigma_true / zeta;
        auto cdf = [](double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); };
        const auto make_rows = [&](double sigma_now) {
            const double s_now = sigma_now / zeta;
            std::vector<MosflmReflection> rows;
            for (int n = 0; n < 50; ++n) {
                const double phic = 10.0 + 6.0 * n + 0.037 * n;
                const double reach = 1.3 / zeta;
                for (int f = static_cast<int>(std::floor((phic - reach) / osc));
                     f * osc <= phic + reach; ++f) {
                    const double mid = f * osc + 0.5 * osc;
                    const double f0 = f * osc - phic;
                    const double f1 = f0 + osc;
                    auto r = slice(f, 0.0f, 0.0f, 0.0f, 0);
                    r.h = r.h0 = n;
                    r.rot = static_cast<float>(phic);
                    r.width = static_cast<float>(6.0 * s_now);
                    r.intensity = static_cast<float>(
                        5000.0 * (cdf(f1 / s_true) - cdf(f0 / s_true))
                    );
                    r.phi_obs = static_cast<float>(mid);
                    rows.push_back(r);
                }
            }
            return rows;
        };
        const auto narrow = fit_rocking_width(make_rows(0.25), 0.25, osc);
        const auto wide = fit_rocking_width(make_rows(0.90), 0.90, osc);
        CHECK(narrow.has_value() && wide.has_value());
        CHECK(narrow.has_value() && close(narrow->sigma_m_deg, sigma_true, 0.05));
        CHECK(wide.has_value() && close(wide->sigma_m_deg, sigma_true, 0.05));
        CHECK(narrow.has_value() && wide.has_value() &&
              std::fabs(narrow->sigma_m_deg - wide->sigma_m_deg) < 0.05);
        const auto few = make_rows(0.25);
        CHECK(!fit_rocking_width({few.begin(), few.begin() + 30}, 0.25, osc).has_value());
    }

    // A constant background on every slice is taken from the wings of the
    // cap and must not widen the second moment.
    {
        using mosflm::diffraction::fit_rocking_width;
        const double osc = 0.2;
        const double zeta = 0.5;
        const double sigma_true = 0.40;
        const double sigma_now = 0.50;
        const double s_true = sigma_true / zeta;
        const double s_now = sigma_now / zeta;
        auto cdf = [](double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); };
        std::vector<MosflmReflection> rows3;
        for (int n = 0; n < 50; ++n) {
            const double phic = 10.0 + 6.0 * n + 0.037 * n;
            const double reach = 1.3 / zeta;
            for (int f = static_cast<int>(std::floor((phic - reach) / osc));
                 f * osc <= phic + reach; ++f) {
                const double f0 = f * osc - phic;
                const double f1 = f0 + osc;
                auto r = slice(f, 0.0f, 0.0f, 0.0f, 0);
                r.h = r.h0 = n;
                r.rot = static_cast<float>(phic);
                r.width = static_cast<float>(6.0 * s_now);
                r.intensity = static_cast<float>(
                    5000.0 * (cdf(f1 / s_true) - cdf(f0 / s_true)) + 20.0
                );
                r.phi_obs = static_cast<float>(f * osc + 0.5 * osc);
                rows3.push_back(r);
            }
        }
        const auto fit = fit_rocking_width(rows3, sigma_now, osc);
        CHECK(fit.has_value() && close(fit->sigma_m_deg, sigma_true, 0.05));
    }

    // Flux beyond |ε| = 1° is outside the sample, so a heavy tail does not
    // move σ_M relative to the same core with those slices removed.
    {
        using mosflm::diffraction::fit_rocking_width;
        const double osc = 0.2;
        const double zeta = 0.5;
        const double sigma_true = 0.40;
        const double sigma_now = 0.50;
        const double s_true = sigma_true / zeta;
        const double s_now = sigma_now / zeta;
        auto cdf = [](double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); };
        const auto lorentz_cdf = [&](double x) {
            return 0.5 + std::atan(x / 1.5) / 3.14159265358979;
        };
        const auto fitted = [&](bool tails) {
            std::vector<MosflmReflection> rows4;
            const double reach = tails ? 3.0 / zeta : 1.0 / zeta;
            for (int n = 0; n < 50; ++n) {
                const double phic = 10.0 + 8.0 * n + 0.037 * n;
                for (int f = static_cast<int>(std::floor((phic - reach) / osc));
                     f * osc <= phic + reach; ++f) {
                    const double f0 = f * osc - phic;
                    const double f1 = f0 + osc;
                    const double eps = zeta * (f0 + f1) * 0.5;
                    auto r = slice(f, 0.0f, 0.0f, 0.0f, 0);
                    r.h = r.h0 = n;
                    r.rot = static_cast<float>(phic);
                    r.width = static_cast<float>(6.0 * s_now);
                    double bin = cdf(f1 / s_true) - cdf(f0 / s_true);
                    if (tails && std::fabs(eps) > 1.0) {
                        bin += 0.15 * (lorentz_cdf(f1) - lorentz_cdf(f0));
                    }
                    r.intensity = static_cast<float>(5000.0 * bin);
                    r.phi_obs = static_cast<float>(f * osc + 0.5 * osc);
                    rows4.push_back(r);
                }
            }
            const auto fit = fit_rocking_width(rows4, sigma_now, osc);
            return fit ? fit->sigma_m_deg : -1.0;
        };
        const double core = fitted(false);
        const double tailed = fitted(true);
        CHECK_MSG(core > 0.0 && tailed > 0.0 && std::fabs(core - tailed) < 0.05,
                  "tails beyond the cap: " + std::to_string(core) + " vs " +
                      std::to_string(tailed));
    }

    // Incomplete edge reflection stays partial.
    const auto edge = sum_partial_slices({slice(1, 0.3f, 30.0f, 0.075f, 0)});
    CHECK(edge.size() == 1u && edge[0].mpart == 1 && close(edge[0].fractioncalc, 0.3, 1e-6));
    // Only a tail slice survived (core frames masked) → reflection dropped.
    const auto tail_only = sum_partial_slices({slice(70, 0.0005f, 0.0f, 10.5f, 0)});
    CHECK(tail_only.empty());
    return mosflm_test::finish("partial_summation");
}
