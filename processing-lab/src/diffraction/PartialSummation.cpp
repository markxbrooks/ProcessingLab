#include "mosflm/diffraction/PartialSummation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

namespace mosflm::diffraction {
namespace {

constexpr int kIbadNegProfile = 4;
/// In-memory only. A neighbour shares the measurement box. Kept so the
/// mosaicity sample can see the peak frame; stripped before the MTZ sum.
constexpr int kFlagSpatialOverlap = 512;

/// True if at least one slice lies in the ±3σ core (contributes intensity).
bool has_core(
    const std::vector<crystal::MosflmReflection>& g,
    const PartialSummationSettings& settings
) {
    return std::any_of(g.begin(), g.end(), [&](const auto& s) {
        return s.fractioncalc >= settings.core_min_fraction &&
               (s.flag & kFlagSpatialOverlap) == 0;
    });
}

crystal::MosflmReflection sum_group(
    const std::vector<crystal::MosflmReflection>& g,
    const PartialSummationSettings& settings
) {
    const auto main = std::max_element(g.begin(), g.end(), [](const auto& a, const auto& b) {
        return a.fractioncalc < b.fractioncalc;
    });
    crystal::MosflmReflection out = *main;
    double i = 0.0, v = 0.0, ipr = 0.0, vpr = 0.0, frac = 0.0, rot = 0.0;
    double pw = 0.0, pphi = 0.0, xw = 0.0, xx = 0.0, xy = 0.0;
    int flag = 0;
    for (const auto& s : g) {
        if ((s.flag & kFlagSpatialOverlap) != 0) {
            continue;
        }
        if (s.fractioncalc >= settings.core_min_fraction) {
            i += s.intensity;
            v += static_cast<double>(s.sigma) * s.sigma;
            ipr += s.ipr;
            vpr += static_cast<double>(s.sigipr) * s.sigipr;
        }
        frac += s.fractioncalc;
        rot += static_cast<double>(s.rot) * s.fractioncalc;
        const double w = std::max(0.0, static_cast<double>(s.intensity));
        if (std::isfinite(s.phi_obs)) {
            pw += w;
            pphi += w * s.phi_obs;
        }
        if (std::isfinite(s.xdet_obs) && std::isfinite(s.ydet_obs)) {
            xw += w;
            xx += w * s.xdet_obs;
            xy += w * s.ydet_obs;
        }
        if (s.fractioncalc >= settings.flag_min_fraction) {
            flag |= s.flag & ~kIbadNegProfile;
        }
    }
    out.intensity = static_cast<float>(i);
    out.sigma = static_cast<float>(std::sqrt(v));
    out.ipr = static_cast<float>(ipr);
    out.sigipr = static_cast<float>(std::sqrt(vpr));
    out.fractioncalc = static_cast<float>(std::min(1.0, frac));
    out.rot = frac > 0.0 ? static_cast<float>(rot / frac) : main->rot;
    out.mpart = frac >= 0.95 ? 0 : 1;
    out.phi_obs = pw > 0.0 ? static_cast<float>(pphi / pw) : main->phi_obs;
    // Refine φ_obs with signed weights in ±2σ about the current estimate
    // (re-centred, ≤3 passes): over a wide core (σ ≈ 1.4° at MOSAIC 1°,
    // ζ 0.3) clipped noise scattered φ_obs 0.3–0.7° and pulled it toward
    // φcalc (76c1b; XDS obs−calc 0.2–0.4°).
    double spacing = std::numeric_limits<double>::infinity();
    for (std::size_t a = 0; a < g.size(); ++a) {
        for (std::size_t b = a + 1; b < g.size(); ++b) {
            const double d = std::fabs(g[a].phi_obs - g[b].phi_obs);
            if (d > 1e-6) {
                spacing = std::min(spacing, d);
            }
        }
    }
    const double half = std::isfinite(spacing)
                            ? std::max(2.0 * static_cast<double>(main->width) / 6.0,
                                       1.5 * spacing)
                            : 0.0;
    for (int pass = 0; pass < 3 && pw > 0.0 && half > 0.0; ++pass) {
        double sw = 0.0;
        double sp = 0.0;
        for (const auto& s : g) {
            if ((s.flag & kFlagSpatialOverlap) != 0) {
                continue;
            }
            if (std::isfinite(s.phi_obs) && std::fabs(s.phi_obs - out.phi_obs) <= half) {
                sw += s.intensity;
                sp += static_cast<double>(s.intensity) * s.phi_obs;
            }
        }
        if (!(sw > 0.0)) {
            break;
        }
        const double next = sp / sw;
        const bool settled = std::fabs(next - out.phi_obs) < 1e-3;
        out.phi_obs = static_cast<float>(next);
        if (settled) {
            break;
        }
    }
    if (pw > 0.0) {
        // Signed (background-subtracted) weights: noise slices in the wide
        // measurement reach average to zero instead of inflating the spread.
        // Clipped weights fed back through MOSAIC → reach and ran away
        // (76c1b: 0.2° → 2° by image 650, then nothing was predicted).
        double m2 = 0.0;
        double sw = 0.0;
        int n_sig = 0;
        for (const auto& s : g) {
            if ((s.flag & kFlagSpatialOverlap) != 0 || !std::isfinite(s.phi_obs)) {
                continue;
            }
            const double w = static_cast<double>(s.intensity);
            const double d = s.phi_obs - out.phi_obs;
            m2 += w * d * d;
            sw += w;
            n_sig += w > 0.0 ? 1 : 0;
        }
        if (n_sig >= 2 && sw > 0.0 && m2 > 0.0) {
            out.phi_obs_sd = static_cast<float>(std::sqrt(m2 / sw));
        }
    }
    if (xw > 0.0) {
        out.xdet_obs = static_cast<float>(xx / xw);
        out.ydet_obs = static_cast<float>(xy / xw);
    }
    // MOSFLM sets this bit only for a profile intensity below −5σ. Slices
    // may each be slightly negative; the test is on the summed IPR.
    out.flag = flag;
    if (out.sigipr > 0.0f && out.ipr < -5.0f * out.sigipr) {
        out.flag |= kIbadNegProfile;
    }
    return out;
}

/// Slices grouped by original HKL and calculated φ (one group per
/// reflection; a gap above ``max_phic_gap_deg`` starts a new one).
std::vector<std::vector<crystal::MosflmReflection>> group_slices(
    std::vector<crystal::MosflmReflection> slices,
    const PartialSummationSettings& settings
) {
    std::sort(slices.begin(), slices.end(), [](const auto& a, const auto& b) {
        return std::tie(a.h0, a.k0, a.l0, a.rot) < std::tie(b.h0, b.k0, b.l0, b.rot);
    });
    std::vector<std::vector<crystal::MosflmReflection>> groups;
    for (auto& s : slices) {
        const bool fresh =
            groups.empty() || s.h0 != groups.back().back().h0 ||
            s.k0 != groups.back().back().k0 || s.l0 != groups.back().back().l0 ||
            std::fabs(s.rot - groups.back().back().rot) > settings.max_phic_gap_deg;
        if (fresh) {
            groups.emplace_back();
        }
        groups.back().push_back(std::move(s));
    }
    return groups;
}

}  // namespace

std::vector<crystal::MosflmReflection> sum_partial_slices(
    std::vector<crystal::MosflmReflection> slices,
    const PartialSummationSettings& settings
) {
    std::vector<crystal::MosflmReflection> out;
    for (const auto& group : group_slices(std::move(slices), settings)) {
        if (has_core(group, settings)) {
            out.push_back(sum_group(group, settings));
        }
    }
    std::stable_sort(out.begin(), out.end(), [](const auto& a, const auto& b) {
        return a.batch < b.batch;
    });
    return out;
}

std::optional<RockingWidthFit> fit_rocking_width(
    const std::vector<crystal::MosflmReflection>& slices,
    double sigma_m_now_deg,
    double osc_deg,
    const PartialSummationSettings& settings,
    std::size_t min_reflections
) {
    if (!(sigma_m_now_deg > 0.0) || !(osc_deg > 0.0)) {
        return std::nullopt;
    }
    // Fixed cap in the Ewald coordinate ε = ζ Δφ. It does not grow with the
    // mosaicity being estimated, so a wide prediction cannot feed itself.
    constexpr double kEpsCapDeg = 1.0;
    std::vector<double> sigmas;
    for (const auto& g : group_slices(slices, settings)) {
        const auto latest = std::max_element(
            g.begin(), g.end(),
            [](const auto& a, const auto& b) { return a.batch < b.batch; }
        );
        // width = 6 σφ and σφ = σ_M / ζ at the prediction that wrote the slice.
        const double sigma_phi = static_cast<double>(latest->width) / 6.0;
        if (!(sigma_phi > 1e-8)) {
            continue;
        }
        const double zeta = sigma_m_now_deg / sigma_phi;
        if (!(zeta > 1e-6)) {
            continue;
        }
        struct Sample {
            double eps;
            double intensity;
            double variance;
        };
        std::vector<Sample> cap;
        cap.reserve(g.size());
        for (const auto& s : g) {
            if (!(s.sigma > 0.0f) || !std::isfinite(s.phi_obs)) {
                continue;
            }
            const double eps = zeta * (static_cast<double>(s.phi_obs) - s.rot);
            if (std::fabs(eps) > kEpsCapDeg) {
                continue;
            }
            cap.push_back(
                {eps, static_cast<double>(s.intensity),
                 static_cast<double>(s.sigma) * s.sigma}
            );
        }
        if (cap.size() < 3) {
            continue;
        }
        const auto moment = [&](double bg) {
            double sum_i = 0.0;
            double sum_ie2 = 0.0;
            double sum_v = 0.0;
            int n_pos = 0;
            for (const auto& s : cap) {
                const double ii = s.intensity - bg;
                if (!(ii > 0.0)) {
                    continue;
                }
                sum_i += ii;
                sum_ie2 += ii * s.eps * s.eps;
                sum_v += s.variance;
                ++n_pos;
            }
            return std::tuple{sum_i, sum_ie2, sum_v, n_pos};
        };
        // First pass with no background, then a pedestal from slices beyond
        // 2.5σ (and at least 0.75°). A fixed |ε|≥0.6° wing still holds the
        // rocking curve when σ_M is a few tenths of a degree, and subtracting
        // it narrows the next estimate.
        const auto raw = moment(0.0);
        const double raw_i = std::get<0>(raw);
        const double raw_v = std::get<2>(raw);
        const int raw_n = std::get<3>(raw);
        if (raw_n < 3 || !(raw_i > 0.0) || !(raw_v > 0.0) ||
            raw_i < 10.0 * std::sqrt(raw_v)) {
            continue;
        }
        const double bin = zeta * osc_deg;
        const double raw_var = std::get<1>(raw) / raw_i - (bin * bin) / 12.0;
        double bg = 0.0;
        if (raw_var > 0.0) {
            const double wing_at =
                std::clamp(2.5 * std::sqrt(raw_var), 0.75, 0.95);
            std::vector<double> wings;
            wings.reserve(cap.size());
            for (const auto& s : cap) {
                if (std::fabs(s.eps) >= wing_at) {
                    wings.push_back(s.intensity);
                }
            }
            if (wings.size() >= 2) {
                const std::size_t mid = wings.size() / 2;
                std::nth_element(
                    wings.begin(),
                    wings.begin() + static_cast<std::ptrdiff_t>(mid),
                    wings.end()
                );
                bg = std::max(0.0, wings[mid]);
            }
        }
        const auto got = moment(bg);
        const double sum_i = std::get<0>(got);
        const double sum_v = std::get<2>(got);
        const int n_pos = std::get<3>(got);
        if (n_pos < 3 || !(sum_i > 0.0) || !(sum_v > 0.0) ||
            sum_i < 10.0 * std::sqrt(sum_v)) {
            continue;
        }
        const double var = std::get<1>(got) / sum_i - (bin * bin) / 12.0;
        if (var > 0.0) {
            sigmas.push_back(std::sqrt(var));
        }
    }
    if (sigmas.size() < min_reflections) {
        return std::nullopt;
    }
    const std::size_t mid = sigmas.size() / 2;
    std::nth_element(
        sigmas.begin(), sigmas.begin() + static_cast<std::ptrdiff_t>(mid),
        sigmas.end()
    );
    return RockingWidthFit{sigmas[mid], sigmas.size()};
}

}  // namespace mosflm::diffraction
