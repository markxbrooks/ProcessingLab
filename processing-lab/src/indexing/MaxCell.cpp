#include "mosflm/indexing/MaxCell.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

namespace mosflm::indexing {
namespace {

/// Root of ``i`` in a union-find forest, halving the path on the way.
std::size_t find_root(std::vector<std::size_t>& parent, std::size_t i) {
    while (parent[i] != i) {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    return i;
}

/// One reflection: the centroid of a chain of spots, and its φ span.
struct Reflection {
    Spot spot;
    double phi_lo{0.0};
    double phi_hi{0.0};
    /// The chain held two spots of one frame: two reflections ran together.
    bool merged_neighbours{false};
};

/// One spot per reflection: chain spots close on the detector and in φ, and
/// take the intensity-weighted centroid of each chain.
std::vector<Reflection> merge_repeats(const SpotList& spots, const MaxCellSettings& settings) {
    const auto& in = spots.spots;
    const double tol_mm = settings.merge_px * spots.header.pixel_mm;
    std::vector<std::size_t> parent(in.size());
    for (std::size_t i = 0; i < in.size(); ++i) {
        parent[i] = i;
    }
    for (std::size_t i = 0; i < in.size(); ++i) {
        for (std::size_t j = i + 1; j < in.size(); ++j) {
            if (std::fabs(in[i].phi_end_deg - in[j].phi_end_deg) <= settings.merge_dphi_deg &&
                std::hypot(in[i].x_mm - in[j].x_mm, in[i].y_mm - in[j].y_mm) <= tol_mm) {
                parent[find_root(parent, i)] = find_root(parent, j);
            }
        }
    }
    struct Sum {
        double w{0.0}, x{0.0}, y{0.0}, phi{0.0}, intensity{0.0};
        double phi_lo{std::numeric_limits<double>::infinity()};
        double phi_hi{-std::numeric_limits<double>::infinity()};
        std::vector<double> frames;
    };
    std::map<std::size_t, Sum> sums;
    for (std::size_t i = 0; i < in.size(); ++i) {
        auto& t = sums[find_root(parent, i)];
        const double w = std::max(in[i].intensity, 1e-6);
        t.w += w;
        t.x += w * in[i].x_mm;
        t.y += w * in[i].y_mm;
        t.phi += w * in[i].phi_end_deg;
        t.intensity += in[i].intensity;
        t.phi_lo = std::min(t.phi_lo, in[i].phi_end_deg);
        t.phi_hi = std::max(t.phi_hi, in[i].phi_end_deg);
        t.frames.push_back(in[i].phi_end_deg);
    }
    std::vector<Reflection> out;
    out.reserve(sums.size());
    for (auto& [root, t] : sums) {
        Reflection r{in[root], t.phi_lo, t.phi_hi, false};
        r.spot.x_mm = t.x / t.w;
        r.spot.y_mm = t.y / t.w;
        r.spot.phi_end_deg = t.phi / t.w;
        r.spot.intensity = t.intensity;
        std::sort(t.frames.begin(), t.frames.end());
        r.merged_neighbours =
            std::adjacent_find(t.frames.begin(), t.frames.end(), [](double a, double b) {
                return std::fabs(a - b) < 1e-6;
            }) != t.frames.end();
        out.push_back(r);
    }
    return out;
}

/// DIALS ``find_overlaps`` stand-in: spots carry no bounding box, so flag
/// chains that swallowed a neighbour, and reflections within ``overlap_px``
/// of another whose φ span overlaps theirs.
std::vector<bool> overlapping(
    const std::vector<Reflection>& refl,
    double pixel_mm,
    const MaxCellSettings& settings
) {
    std::vector<bool> flag(refl.size(), false);
    const double tol_mm = settings.overlap_px * pixel_mm;
    const double pad = settings.merge_dphi_deg;
    for (std::size_t i = 0; i < refl.size(); ++i) {
        flag[i] = flag[i] || refl[i].merged_neighbours;
        for (std::size_t j = i + 1; j < refl.size(); ++j) {
            const bool phi_overlap = refl[i].phi_lo - pad <= refl[j].phi_hi &&
                                     refl[j].phi_lo - pad <= refl[i].phi_hi;
            if (phi_overlap && std::hypot(refl[i].spot.x_mm - refl[j].spot.x_mm,
                                          refl[i].spot.y_mm - refl[j].spot.y_mm) <= tol_mm) {
                flag[i] = true;
                flag[j] = true;
            }
        }
    }
    return flag;
}

}  // namespace

std::vector<double> ice_ring_d_star_sq(double d_min_A) {
    // Hexagonal metric: 1/d² = 4/(3a²)(h² + hk + k²) + l²/c².
    constexpr double a = 4.498;
    constexpr double c = 7.338;
    std::vector<double> rings;
    if (!(d_min_A > 0.0)) {
        return rings;
    }
    const double s2_max = 1.0 / (d_min_A * d_min_A);
    const int hmax = static_cast<int>(std::ceil(a / d_min_A)) + 1;
    const int lmax = static_cast<int>(std::ceil(c / d_min_A));
    for (int h = -hmax; h <= hmax; ++h) {
        for (int k = -hmax; k <= hmax; ++k) {
            for (int l = 0; l <= lmax; ++l) {
                if (h == 0 && k == 0 && l == 0) {
                    continue;
                }
                // P6₃/mmc: hh(-2h)l and its permutations (with 000l) need l even.
                if (l % 2 != 0 && (h == k || k == -2 * h || h == -2 * k)) {
                    continue;
                }
                const double s2 = 4.0 / (3.0 * a * a) * (h * h + h * k + k * k) + l * l / (c * c);
                if (s2 <= s2_max) {
                    rings.push_back(s2);
                }
            }
        }
    }
    std::sort(rings.begin(), rings.end());
    rings.erase(std::unique(rings.begin(), rings.end(),
                            [](double p, double q) { return std::fabs(p - q) < 1e-9; }),
                rings.end());
    return rings;
}

MaxCellEstimate estimate_max_cell(
    const std::vector<float>& x,
    const std::vector<float>& y,
    const std::vector<float>& z,
    const std::vector<int>& group,
    const MaxCellSettings& settings
) {
    MaxCellEstimate out;
    const std::size_t n = std::min({x.size(), y.size(), z.size(), group.size()});
    std::map<int, std::vector<std::size_t>> members;
    for (std::size_t i = 0; i < n; ++i) {
        members[group[i]].push_back(i);
    }

    // Nearest neighbour of each point within its group, as a direct distance.
    std::vector<double> direct;
    direct.reserve(n);
    for (const auto& [g, idx] : members) {
        for (std::size_t a : idx) {
            double best = std::numeric_limits<double>::infinity();
            for (std::size_t b : idx) {
                if (a == b) {
                    continue;
                }
                const double dx = x[a] - x[b];
                const double dy = y[a] - y[b];
                const double dz = z[a] - z[b];
                best = std::min(best, dx * dx + dy * dy + dz * dz);
            }
            if (best > 0.0 && std::isfinite(best)) {
                direct.push_back(1.0 / std::sqrt(best));
            }
        }
    }
    // Drop nonsensical distances; the histogram rule also drops the top 1%.
    std::erase_if(direct, [](double d) { return !(d > 1.0); });
    std::sort(direct.begin(), direct.end());
    constexpr std::size_t kNear = 10;
    if (direct.size() <= kNear) {
        return out;
    }
    if (settings.percentile > 0.0) {
        const auto from_top = static_cast<std::size_t>((1.0 - std::min(settings.percentile, 1.0)) *
                                                       static_cast<double>(direct.size()));
        out.ok = true;
        out.max_cell_A = settings.tolerance * direct[direct.size() - 1 - std::min(from_top, direct.size() - 1)];
        out.n_neighbours = direct.size();
        return out;
    }
    direct.resize(static_cast<std::size_t>(std::floor(0.99 * static_cast<double>(direct.size()))));
    const int n_slots = static_cast<int>(direct.size()) / std::max(1, settings.nn_per_bin);
    if (direct.size() <= kNear || n_slots < 1) {
        return out;
    }

    // Histogram over the data range (linear or log10). Heights are counts per
    // Å so log bins compare fairly; the last bin above the height fraction of
    // the tallest sets the most probable large separation.
    const auto to_axis = [&](double d) { return settings.log_binning ? std::log10(d) : d; };
    const auto from_axis = [&](double u) { return settings.log_binning ? std::pow(10.0, u) : u; };
    const double lo = to_axis(direct.front());
    const double width = std::max(to_axis(direct.back()) - lo, 1e-9) / n_slots;
    std::vector<int> counts(static_cast<std::size_t>(n_slots), 0);
    for (double d : direct) {
        const int k = std::min(n_slots - 1, static_cast<int>((to_axis(d) - lo) / width));
        ++counts[static_cast<std::size_t>(k)];
    }
    std::vector<double> height(counts.size());
    for (int k = 0; k < n_slots; ++k) {
        const double span = from_axis(lo + (k + 1) * width) - from_axis(lo + k * width);
        height[static_cast<std::size_t>(k)] = counts[static_cast<std::size_t>(k)] / span;
    }
    const double highest = *std::max_element(height.begin(), height.end());
    int last = 0;
    for (int k = 0; k < n_slots; ++k) {
        if (height[static_cast<std::size_t>(k)] > settings.max_height_fraction * highest) {
            last = k;
        }
    }
    out.ok = true;
    out.max_cell_A = settings.tolerance * from_axis(lo + (last + 1) * width);
    out.n_neighbours = direct.size();
    return out;
}

MaxCellEstimate estimate_max_cell(
    const SpotList& spots,
    const SpotToRlpGeometry& geometry,
    const MaxCellSettings& settings
) {
    std::vector<Reflection> refl = merge_repeats(spots, settings);
    const std::size_t n_reflections = refl.size();
    SpotList merged;
    merged.header = spots.header;
    merged.spots.reserve(refl.size());
    for (const auto& r : refl) {
        merged.spots.push_back(r.spot);
    }
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> z;
    if (spots_to_rlp_xyz(merged, {}, geometry, x, y, z) != merged.spots.size()) {
        return {};
    }

    // DIALS rejects ice-ring and overlapping reflections before the analysis.
    std::vector<bool> reject(refl.size(), false);
    MaxCellEstimate out;
    if (settings.filter_ice && !refl.empty()) {
        double s2_max = 0.0;
        std::vector<double> s2(refl.size());
        for (std::size_t i = 0; i < refl.size(); ++i) {
            s2[i] = static_cast<double>(x[i]) * x[i] + static_cast<double>(y[i]) * y[i] +
                    static_cast<double>(z[i]) * z[i];
            s2_max = std::max(s2_max, s2[i]);
        }
        const double half = 0.5 * settings.ice_width;
        const auto rings = ice_ring_d_star_sq(1.0 / std::sqrt(s2_max + half));
        for (std::size_t i = 0; i < refl.size(); ++i) {
            const auto it = std::lower_bound(rings.begin(), rings.end(), s2[i] - half);
            if (it != rings.end() && *it < s2[i] + half) {
                reject[i] = true;
                ++out.n_ice;
            }
        }
    }
    if (settings.filter_overlaps) {
        const auto flag = overlapping(refl, spots.header.pixel_mm, settings);
        for (std::size_t i = 0; i < refl.size(); ++i) {
            if (flag[i] && !reject[i]) {
                reject[i] = true;
                ++out.n_overlaps;
            }
        }
    }

    // Group by φ block and by side of the rotation axis (y in the spot→RLP
    // frame), as DIALS separates entering from exiting reflections.
    constexpr double kDeg2Rad = 3.14159265358979323846 / 180.0;
    const double omega = spots.header.omega_deg * kDeg2Rad;
    const double finvertx = (spots.header.invflag != 0) ? -1.0 : 1.0;
    double phi_min = std::numeric_limits<double>::infinity();
    for (const auto& s : merged.spots) {
        phi_min = std::min(phi_min, s.phi_end_deg);
    }
    const double step = settings.step_deg > 0.0 ? settings.step_deg : 360.0;
    std::vector<float> kx;
    std::vector<float> ky;
    std::vector<float> kz;
    std::vector<int> group;
    for (std::size_t i = 0; i < merged.spots.size(); ++i) {
        if (reject[i]) {
            continue;
        }
        const auto& s = merged.spots[i];
        const double xt = s.x_mm - geometry.beam_x_mm;
        const double yt = s.y_mm - geometry.beam_y_mm;
        const double xs = xt * std::cos(omega) * finvertx + yt * std::sin(omega);
        const int block = static_cast<int>(std::floor((s.phi_end_deg - phi_min) / step));
        kx.push_back(x[i]);
        ky.push_back(y[i]);
        kz.push_back(z[i]);
        group.push_back(2 * block + (xs >= 0.0 ? 1 : 0));
    }
    const MaxCellEstimate nn = estimate_max_cell(kx, ky, kz, group, settings);
    out.ok = nn.ok;
    out.max_cell_A = nn.max_cell_A;
    out.n_neighbours = nn.n_neighbours;
    out.n_reflections = n_reflections;
    return out;
}

}  // namespace mosflm::indexing
