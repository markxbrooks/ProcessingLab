#include "mosflm/diffraction/PartialSummation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

namespace mosflm::diffraction {
namespace {

constexpr int kIbadNegProfile = 4;

/// True if at least one slice lies in the ±3σ core (contributes intensity).
bool has_core(
    const std::vector<crystal::MosflmReflection>& g,
    const PartialSummationSettings& settings
) {
    return std::any_of(g.begin(), g.end(), [&](const auto& s) {
        return s.fractioncalc >= settings.core_min_fraction;
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
            if (!std::isfinite(s.phi_obs)) {
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
    out.flag = flag | (out.ipr < 0.0f ? kIbadNegProfile : 0);
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
    double width_now_deg,
    double osc_deg,
    const PartialSummationSettings& settings,
    std::size_t min_reflections
) {
    if (!(width_now_deg > 0.0) || !(osc_deg > 0.0)) {
        return std::nullopt;
    }
    // Per slice: frame φ range relative to φcalc, intensity, 1/σ².
    struct Slice {
        double lo, hi, i, w;
    };
    struct Refl {
        double sigma_now;  // σφ at width_now_deg
        std::vector<Slice> s;
    };
    std::vector<Refl> refl;
    for (const auto& g : group_slices(slices, settings)) {
        if (!has_core(g, settings)) {
            continue;
        }
        double isum = 0.0, vsum = 0.0, frac = 0.0;
        // The window spans earlier periods predicted with another MOSAIC;
        // the latest slice was predicted with width_now_deg.
        const auto latest = std::max_element(g.begin(), g.end(), [](const auto& a, const auto& b) {
            return a.batch < b.batch;
        });
        Refl r{static_cast<double>(latest->width) / 6.0, {}};
        double mw = 0.0;
        for (const auto& s : g) {
            if (!(s.sigma > 0.0f) || !std::isfinite(s.phi_obs) || s.flag != 0) {
                continue;
            }
            const double mid = static_cast<double>(s.phi_obs) - s.rot;
            r.s.push_back({mid - 0.5 * osc_deg, mid + 0.5 * osc_deg, s.intensity,
                           1.0 / (static_cast<double>(s.sigma) * s.sigma)});
            isum += s.intensity;
            vsum += static_cast<double>(s.sigma) * s.sigma;
            frac += s.fractioncalc;
            mw += s.intensity * mid;
        }
        // Centre the model on the observed φ centroid: an orientation error
        // is for positional refinement, not extra rocking width.
        if (isum > 0.0) {
            const double shift = mw / isum;
            for (auto& sl : r.s) {
                sl.lo -= shift;
                sl.hi -= shift;
            }
        }
        // Strong and (under the current model) fully measured.
        if (r.sigma_now > 0.0 && r.s.size() >= 2 && frac >= 0.95 &&
            isum >= 10.0 * std::sqrt(vsum)) {
            refl.push_back(std::move(r));
        }
    }
    if (refl.size() < min_reflections) {
        return std::nullopt;
    }
    auto cdf = [](double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); };
    // Σ over reflections of the unexplained fraction of Σ w I² after fitting
    // the model profile with a free amplitude and a free per-slice baseline.
    // The baseline absorbs a background bias in the per-frame intensities:
    // without it a slight positive bias in the outer slices reads as heavier
    // tails, the fit widens, the wider model adds more outer slices, and the
    // width runs away (76c1b: 0.96° → 2.0°, XDS 1.04°).
    auto cost = [&](double width) {
        double c = 0.0;
        for (const auto& r : refl) {
            const double sig = r.sigma_now * width / width_now_deg;
            double sw = 0.0, sp = 0.0, si = 0.0, spp = 0.0, sip = 0.0, sii = 0.0;
            for (const auto& s : r.s) {
                const double pj = cdf(s.hi / sig) - cdf(s.lo / sig);
                sw += s.w;
                sp += s.w * pj;
                si += s.w * s.i;
                spp += s.w * pj * pj;
                sip += s.w * s.i * pj;
                sii += s.w * s.i * s.i;
            }
            if (!(sii > 0.0) || !(sw > 0.0)) {
                continue;
            }
            // Weighted least squares for I ≈ t·p + b; residual Σ w (I − t p − b)².
            const double cpp = spp - sp * sp / sw;
            const double cip = sip - si * sp / sw;
            const double cii = sii - si * si / sw;
            const double resid = cpp > 0.0 ? cii - cip * cip / cpp : cii;
            c += std::max(resid, 0.0) / sii;
        }
        return c;
    };
    // Coarse log grid, then golden-section between the best's neighbours.
    constexpr int kGrid = 48;
    const double x0 = std::log(0.02);
    const double x1 = std::log(3.0);
    const double dx = (x1 - x0) / (kGrid - 1);
    int best = 0;
    double best_c = std::numeric_limits<double>::infinity();
    for (int k = 0; k < kGrid; ++k) {
        const double c = cost(std::exp(x0 + dx * k));
        if (c < best_c) {
            best_c = c;
            best = k;
        }
    }
    double a = x0 + dx * std::max(best - 1, 0);
    double b = x0 + dx * std::min(best + 1, kGrid - 1);
    const double gr = 0.5 * (std::sqrt(5.0) - 1.0);
    double c1 = b - gr * (b - a);
    double c2 = a + gr * (b - a);
    double f1 = cost(std::exp(c1));
    double f2 = cost(std::exp(c2));
    for (int it = 0; it < 40 && b - a > 1e-4; ++it) {
        if (f1 < f2) {
            b = c2;
            c2 = c1;
            f2 = f1;
            c1 = b - gr * (b - a);
            f1 = cost(std::exp(c1));
        } else {
            a = c1;
            c1 = c2;
            f1 = f2;
            c2 = a + gr * (b - a);
            f2 = cost(std::exp(c2));
        }
    }
    return RockingWidthFit{std::exp(0.5 * (a + b)), refl.size()};
}

}  // namespace mosflm::diffraction
