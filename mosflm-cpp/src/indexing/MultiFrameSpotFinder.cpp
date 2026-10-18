#include "mosflm/indexing/MultiFrameSpotFinder.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace mosflm::indexing {
namespace {

struct Acc {
    double x_sum{0.0};
    double y_sum{0.0};
    double phi_sum{0.0};
    double inten_sum{0.0};
    int n{0};
};

double mid_phi(const Spot& s) {
    return 0.5 * (s.phi_start_deg + s.phi_end_deg);
}

}  // namespace

SpotList find_spots_multiframe(
    int first_frame,
    int last_frame,
    const FrameLoader& load,
    const MultiFrameSpotSettings& settings
) {
    SpotFinder finder;
    std::vector<Spot> raw;
    SpotListHeader header;
    for (int f = first_frame; f <= last_frame; ++f) {
        const auto image = load(f);
        SpotFindSettings fs = settings.per_frame;
        if (fs.phi_end_deg <= fs.phi_start_deg &&
            image.metadata.oscillation_range_deg > 0.0) {
            fs.phi_start_deg = image.metadata.oscillation_start_deg;
            fs.phi_end_deg =
                fs.phi_start_deg + image.metadata.oscillation_range_deg;
        }
        auto list = finder.find(image, fs);
        if (f == first_frame) {
            header = list.header;
        }
        for (auto& s : list.spots) {
            raw.push_back(std::move(s));
        }
    }

    std::sort(
        raw.begin(), raw.end(),
        [](const Spot& a, const Spot& b) { return a.intensity > b.intensity; }
    );

    std::vector<Acc> clusters;
    std::vector<char> used(raw.size(), 0);
    const double xy2 = settings.merge_xy_mm * settings.merge_xy_mm;
    for (std::size_t i = 0; i < raw.size(); ++i) {
        if (used[i]) {
            continue;
        }
        Acc acc;
        const auto& seed = raw[i];
        const double seed_phi = mid_phi(seed);
        for (std::size_t j = i; j < raw.size(); ++j) {
            if (used[j]) {
                continue;
            }
            const auto& s = raw[j];
            const double dx = s.x_mm - seed.x_mm;
            const double dy = s.y_mm - seed.y_mm;
            const double dphi = std::fabs(mid_phi(s) - seed_phi);
            if (dx * dx + dy * dy > xy2) {
                continue;
            }
            if (dphi > settings.merge_phi_tol_deg) {
                continue;
            }
            used[j] = 1;
            acc.x_sum += s.x_mm;
            acc.y_sum += s.y_mm;
            acc.phi_sum += mid_phi(s);
            acc.inten_sum += s.intensity;
            ++acc.n;
        }
        if (acc.n >= settings.min_frames) {
            clusters.push_back(acc);
        }
    }

    std::sort(
        clusters.begin(), clusters.end(),
        [](const Acc& a, const Acc& b) { return a.inten_sum > b.inten_sum; }
    );
    if (static_cast<int>(clusters.size()) > settings.max_spots) {
        clusters.resize(static_cast<std::size_t>(settings.max_spots));
    }

    SpotList out;
    out.header = header;
    out.spots.reserve(clusters.size());
    for (const auto& c : clusters) {
        Spot s;
        const double inv = 1.0 / static_cast<double>(c.n);
        s.x_mm = c.x_sum * inv;
        s.y_mm = c.y_sum * inv;
        const double phi = c.phi_sum * inv;
        s.phi_start_deg = phi;
        s.phi_end_deg = phi;
        s.intensity = c.inten_sum;
        out.spots.push_back(s);
    }
    return out;
}

}  // namespace mosflm::indexing
