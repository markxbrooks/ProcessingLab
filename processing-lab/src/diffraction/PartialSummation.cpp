#include "mosflm/diffraction/PartialSummation.hpp"

#include <algorithm>
#include <cmath>
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
    if (pw > 0.0) {
        double m2 = 0.0;
        int n_sig = 0;
        for (const auto& s : g) {
            const double w = std::max(0.0, static_cast<double>(s.intensity));
            if (w > 0.0 && std::isfinite(s.phi_obs)) {
                const double d = s.phi_obs - out.phi_obs;
                m2 += w * d * d;
                ++n_sig;
            }
        }
        if (n_sig >= 2) {
            out.phi_obs_sd = static_cast<float>(std::sqrt(m2 / pw));
        }
    }
    if (xw > 0.0) {
        out.xdet_obs = static_cast<float>(xx / xw);
        out.ydet_obs = static_cast<float>(xy / xw);
    }
    out.flag = flag | (out.ipr < 0.0f ? kIbadNegProfile : 0);
    return out;
}

}  // namespace

std::vector<crystal::MosflmReflection> sum_partial_slices(
    std::vector<crystal::MosflmReflection> slices,
    const PartialSummationSettings& settings
) {
    std::sort(slices.begin(), slices.end(), [](const auto& a, const auto& b) {
        return std::tie(a.h0, a.k0, a.l0, a.rot) < std::tie(b.h0, b.k0, b.l0, b.rot);
    });
    std::vector<crystal::MosflmReflection> out;
    std::vector<crystal::MosflmReflection> group;
    for (auto& s : slices) {
        if (!group.empty()) {
            const auto& last = group.back();
            if (s.h0 != last.h0 || s.k0 != last.k0 || s.l0 != last.l0 ||
                std::fabs(s.rot - last.rot) > settings.max_phic_gap_deg) {
                if (has_core(group, settings)) {
                    out.push_back(sum_group(group, settings));
                }
                group.clear();
            }
        }
        group.push_back(std::move(s));
    }
    if (!group.empty() && has_core(group, settings)) {
        out.push_back(sum_group(group, settings));
    }
    std::stable_sort(out.begin(), out.end(), [](const auto& a, const auto& b) {
        return a.batch < b.batch;
    });
    return out;
}

}  // namespace mosflm::diffraction
