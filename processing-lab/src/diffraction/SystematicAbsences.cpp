#include "mosflm/diffraction/SystematicAbsences.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <tuple>

namespace mosflm::diffraction {
namespace {

constexpr double kPi = 3.14159265358979323846;
using Hkl = std::tuple<int, int, int>;

Hkl p4_asu(int h, int k, int l) {
    std::array<Hkl, 8> eq{{{h,k,l}, {-k,h,l}, {-h,-k,l}, {k,-h,l},
                           {-h,-k,-l}, {k,-h,-l}, {h,k,-l}, {-k,h,-l}}};
    return *std::min_element(eq.begin(), eq.end());
}

std::pair<double, int> twofold_cc(
    const std::map<Hkl, double>& intensities, bool about_a
) {
    std::set<std::pair<Hkl, Hkl>> seen;
    std::vector<double> x, y;
    for (const auto& [hkl, intensity] : intensities) {
        const auto [h, k, l] = hkl;
        const Hkl other = about_a ? p4_asu(h, -k, -l)
                                  : p4_asu(-h, k, -l);
        if (other == hkl || !intensities.contains(other))
            continue;
        auto key = std::minmax(hkl, other);
        if (!seen.emplace(key.first, key.second).second)
            continue;
        x.push_back(intensity);
        y.push_back(intensities.at(other));
    }
    if (x.size() < 3)
        return {0.0, static_cast<int>(x.size())};
    const double mx = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    const double my = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
    double cov = 0.0, vx = 0.0, vy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double dx = x[i] - mx, dy = y[i] - my;
        cov += dx * dy; vx += dx * dx; vy += dy * dy;
    }
    return {vx > 0.0 && vy > 0.0 ? cov / std::sqrt(vx * vy) : 0.0,
            static_cast<int>(x.size())};
}

AbsenceZoneScore score_axis(
    const std::vector<crystal::MosflmReflection>& reflections,
    int axis,
    bool pool_ab = false
) {
    std::map<int, std::pair<double, int>> values;
    for (const auto& r : reflections) {
        const bool on_a = r.k == 0 && r.l == 0;
        const bool on_b = r.h == 0 && r.l == 0;
        const bool on_c = r.h == 0 && r.k == 0;
        int index = 0;
        if (pool_ab && axis != 2) {
            // Tetragonal h00 and 0k0 are one zone; the P4 ASU stores it on 0k0.
            if (on_a)
                index = r.h;
            else if (on_b)
                index = r.k;
            else
                continue;
        } else if (axis == 0 && on_a) {
            index = r.h;
        } else if (axis == 1 && on_b) {
            index = r.k;
        } else if (axis == 2 && on_c) {
            index = r.l;
        } else {
            continue;
        }
        if (index == 0 || r.sigma <= 0.0f || !std::isfinite(r.intensity) ||
            !std::isfinite(r.sigma))
            continue;
        auto& v = values[std::abs(index)];
        v.first += static_cast<double>(r.intensity / r.sigma);
        ++v.second;
    }
    for (auto& [i, v] : values) {
        (void)i;
        v.first /= static_cast<double>(v.second);
    }

    AbsenceZoneScore best;
    best.axis = "abc"[axis];
    best.observations = static_cast<int>(values.size());
    if (values.size() < 6)
        return best;

    // POINTLESS-style I' adjustment reduces contamination from a nearby
    // strong axial reflection before evaluating periodicity.
    std::map<int, double> adjusted;
    for (const auto& [i, v] : values) {
        const double left = values.contains(i - 1) ? values.at(i - 1).first : 0.0;
        const double right = values.contains(i + 1) ? values.at(i + 1).first : 0.0;
        adjusted[i] = v.first - 0.02 * (left + right);
    }
    for (int period : {2, 3, 4, 6}) {
        // Fundamental harmonic distinguishes a 4_1/4_3 pattern (l=4n)
        // from 4_2 (l=2n), whose signal is at the second harmonic.
        for (int harmonic : {1}) {
            if (harmonic >= period)
                continue;
            double re = 0.0, im = 0.0, origin = 0.0;
            for (const auto& [i, x] : adjusted) {
                const double phi = 2.0 * kPi * harmonic * i / period;
                re += x * std::cos(phi);
                im += x * std::sin(phi);
                origin += std::abs(x);
            }
            if (origin <= 0.0)
                continue;
            const double peak = std::hypot(re, im) / origin;
            if (peak < 1e-9)
                continue;
            const double sd = 1.0 / std::sqrt(2.0 * static_cast<double>(values.size()));
            const double z = peak / sd;
            const double probability = std::erf(z / std::sqrt(2.0));
            const double best_z = best.peak_height / std::max(best.sd, 1e-12);
            if (z > best_z + 1e-9 ||
                (std::fabs(z - best_z) <= 1e-9 && period > best.period)) {
                best.period = period;
                best.peak_height = peak;
                best.sd = sd;
                best.probability = probability;
            }
        }
    }
    return best;
}

}  // namespace

AbsenceAnalysis analyze_systematic_absences(
    const std::vector<crystal::MosflmReflection>& reflections,
    bool tetragonal_metric
) {
    AbsenceAnalysis out;
    for (int i = 0; i < 3; ++i)
        out.axes[static_cast<std::size_t>(i)] =
            score_axis(reflections, i, tetragonal_metric);

    // A tetragonal axial signature matching 4_1(c) and 2_1(a,b), plus
    // intensity agreement under both added twofold operators, supports the
    // P 41 21 2 / P 43 21 2 pair.
    const auto strong = [](const AbsenceZoneScore& z, int period) {
        return z.period == period && z.observations >= 8 && z.probability >= 0.95;
    };
    if (tetragonal_metric && strong(out.axes[0], 2) &&
        strong(out.axes[1], 2) && strong(out.axes[2], 4)) {
        std::map<Hkl, std::pair<double, int>> sums;
        for (const auto& r : reflections) {
            if (r.sigma <= 0.0f || !std::isfinite(r.intensity) ||
                !std::isfinite(r.sigma) || (r.h == 0 && r.k == 0 && r.l == 0))
                continue;
            auto& v = sums[p4_asu(r.h, r.k, r.l)];
            v.first += r.intensity;
            ++v.second;
        }
        std::map<Hkl, double> means;
        for (const auto& [hkl, v] : sums)
            means[hkl] = v.first / v.second;
        const auto a = twofold_cc(means, true);
        const auto b = twofold_cc(means, false);
        out.twofold_a_cc = a.first;
        out.twofold_a_pairs = a.second;
        out.twofold_b_cc = b.first;
        out.twofold_b_pairs = b.second;
        out.laue_group_checked = a.second >= 20 && b.second >= 20;
        if (out.laue_group_checked && a.first >= 0.6 && b.first >= 0.6) {
            out.suggested_space_group = 92;  // #96 is the indistinguishable enantiomorph
            out.suggestion_confidence = std::min({out.axes[0].probability,
                                                  out.axes[1].probability,
                                                  out.axes[2].probability,
                                                  a.first, b.first});
        }
    }
    return out;
}

}  // namespace mosflm::diffraction
