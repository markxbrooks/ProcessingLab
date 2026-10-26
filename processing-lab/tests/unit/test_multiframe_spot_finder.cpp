/**
 * MULTIFRAME spot merging: tracks are one-to-one, unambiguous and anchored
 * to their centroid; settings are validated; ordering is deterministic.
 */

#include "mosflm/indexing/MultiFrameSpotFinder.hpp"
#include "test_util.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

using mosflm::image::DiffractionImage;
using mosflm::indexing::find_spots_multiframe;
using mosflm::indexing::MultiFrameSpotSettings;
using mosflm::indexing::SpotList;

namespace {

constexpr int kSize = 64;
constexpr double kPixelMm = 0.1;

/// Bright pixels per frame: frame → [(col, row, counts)].
using Scene = std::map<int, std::vector<std::tuple<int, int, int>>>;

DiffractionImage frame_of(const Scene& scene, int f) {
    DiffractionImage img;
    img.metadata.width = kSize;
    img.metadata.height = kSize;
    img.metadata.pixel_size_x_mm = kPixelMm;
    img.metadata.pixel_size_y_mm = kPixelMm;
    img.metadata.beam_x_mm = 3.2;
    img.metadata.beam_y_mm = 3.2;
    img.metadata.oscillation_start_deg = 0.2 * (f - 1);
    img.metadata.oscillation_range_deg = 0.2;
    img.pixels.assign(static_cast<std::size_t>(kSize * kSize), 10);
    if (const auto it = scene.find(f); it != scene.end()) {
        for (const auto& [c, r, v] : it->second) {
            img.pixels[static_cast<std::size_t>(r) * kSize + static_cast<std::size_t>(c)] = v;
        }
    }
    return img;
}

MultiFrameSpotSettings settings(double radius_px) {
    MultiFrameSpotSettings s;
    s.per_frame.threshold_sigma = 5.0f;
    s.per_frame.min_distance = 4;
    s.per_frame.max_spots = 50;
    s.per_frame.border = 2;
    s.merge_xy_mm = radius_px * kPixelMm;
    s.max_spots = 100;
    return s;
}

SpotList run(const Scene& scene, int last, const MultiFrameSpotSettings& s) {
    return find_spots_multiframe(1, last, [&](int f) { return frame_of(scene, f); }, s);
}

bool throws_invalid(int first, int last, const MultiFrameSpotSettings& s) {
    try {
        (void)find_spots_multiframe(first, last, [](int) { return DiffractionImage{}; }, s);
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

}  // namespace

int main() {
    // One reflection on five frames with ±1 px jitter → one spot.
    {
        Scene scene{{1, {{20, 20, 800}}}, {2, {{21, 20, 2000}}}, {3, {{20, 21, 3000}}},
                    {4, {{20, 20, 2000}}}, {5, {{19, 20, 800}}}};
        const auto out = run(scene, 5, settings(2.5));
        CHECK_EQ(out.spots.size(), std::size_t{1});
        CHECK_CLOSE(out.spots[0].intensity, 3000.0, 1.0);  // the brightest frame's peak
    }

    // Two spots on one frame both within reach of one spot on the next: the
    // match is ambiguous, so nothing is joined (three spots, not one).
    {
        Scene scene{{1, {{20, 30, 2000}, {28, 30, 2000}}}, {2, {{24, 30, 3000}}}};
        const auto out = run(scene, 2, settings(5.0));
        CHECK_EQ(out.spots.size(), std::size_t{3});
    }

    // Two reflections side by side on several frames stay two tracks.
    {
        Scene scene;
        for (int f = 1; f <= 4; ++f) {
            scene[f] = {{20, 30, 1000 + 500 * f}, {28, 30, 1200 + 400 * f}};
        }
        const auto out = run(scene, 4, settings(2.5));
        CHECK_EQ(out.spots.size(), std::size_t{2});
    }

    // A chain drifting 2 px per frame: each step is within the 2.5 px radius,
    // but the track is anchored to its centroid (4 radii = 10 px), so a
    // 26 px run is cut rather than followed to the end.
    {
        Scene scene;
        for (int f = 1; f <= 14; ++f) {
            scene[f] = {{6 + 2 * (f - 1), 30, 2000}};
        }
        const auto out = run(scene, 14, settings(2.5));
        CHECK(out.spots.size() >= 2);
        // A short drift (3 px over 4 frames) is one reflection, followed.
        Scene short_drift;
        for (int f = 1; f <= 4; ++f) {
            short_drift[f] = {{20 + (f - 1), 30, 2000}};
        }
        CHECK_EQ(run(short_drift, 4, settings(2.5)).spots.size(), std::size_t{1});
    }

    // A frame without the reflection ends its track: two spots.
    {
        Scene scene{{1, {{20, 20, 2000}}}, {2, {{20, 20, 2000}}}, {4, {{20, 20, 2000}}}, {5, {{20, 20, 2000}}}};
        const auto out = run(scene, 5, settings(2.5));
        CHECK_EQ(out.spots.size(), std::size_t{2});
        // min_frames 3 drops both two-frame tracks.
        auto s = settings(2.5);
        s.min_frames = 3;
        CHECK(run(scene, 5, s).spots.empty());
    }

    // max_spots ≤ 0 is unlimited in both paths (it used to resize the merged
    // list to a huge unsigned size).
    {
        Scene scene{{1, {{10, 10, 2000}, {30, 30, 3000}, {50, 50, 4000}}}, {2, {{10, 50, 2500}}}};
        for (const int max_spots : {0, -1}) {
            auto s = settings(2.5);
            s.max_spots = max_spots;
            CHECK_EQ(run(scene, 2, s).spots.size(), std::size_t{4});
            s.merge_xy_mm = 0.0;
            CHECK_EQ(run(scene, 2, s).spots.size(), std::size_t{4});
        }
        auto s = settings(2.5);
        s.max_spots = 2;
        const auto two = run(scene, 2, s);
        CHECK_EQ(two.spots.size(), std::size_t{2});
        CHECK_CLOSE(two.spots[0].intensity, 4000.0, 1.0);  // brightest first
    }

    // Invalid input is rejected.
    {
        CHECK(throws_invalid(5, 4, settings(2.5)));
        auto s = settings(2.5);
        s.merge_xy_mm = std::numeric_limits<double>::quiet_NaN();
        CHECK(throws_invalid(1, 2, s));
        s.merge_xy_mm = -0.1;
        CHECK(throws_invalid(1, 2, s));
        s = settings(2.5);
        s.min_frames = 0;
        CHECK(throws_invalid(1, 2, s));
    }

    // Serial and parallel give the same list, in the same order, even with
    // equal intensities.
    {
        Scene scene;
        for (int f = 1; f <= 12; ++f) {
            scene[f] = {{8 + 4 * (f % 5), 10, 2000}, {40, 8 + 4 * (f % 7), 2000}, {30, 50, 1500 + 100 * (f % 3)}};
        }
        auto serial = settings(2.5);
        auto parallel = serial;
        parallel.nproc = 4;
        const auto a = run(scene, 12, serial);
        const auto b = run(scene, 12, parallel);
        CHECK_EQ(a.spots.size(), b.spots.size());
        for (std::size_t i = 0; i < std::min(a.spots.size(), b.spots.size()); ++i) {
            CHECK(a.spots[i].x_mm == b.spots[i].x_mm && a.spots[i].y_mm == b.spots[i].y_mm &&
                  a.spots[i].phi_end_deg == b.spots[i].phi_end_deg);
        }
    }

    return mosflm_test::finish("multiframe_spot_finder");
}
