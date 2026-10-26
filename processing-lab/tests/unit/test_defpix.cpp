#include "mosflm/control/ControlSession.hpp"
#include "mosflm/diffraction/Defpix.hpp"
#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/image/DiffractionImage.hpp"
#include "mosflm/io/XdsInpImporter.hpp"
#include "test_util.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using mosflm::diffraction::PixelStatus;
using mosflm::image::DiffractionImage;
using mosflm::image::kPixelMasked;
using mosflm::image::kPixelOverloaded;

namespace {

constexpr int kW = 200;
constexpr int kH = 200;

std::size_t at(int c, int r) {
    return static_cast<std::size_t>(r) * kW + static_cast<std::size_t>(c);
}

/// Frame ``f`` (0–9) of a synthetic sweep: background ~20 counts, a module
/// gap (row 0), a stuck pixel, a flickering one, a hot one, a spot that
/// saturates one pixel in one frame only and a shadow at 10 % background.
DiffractionImage frame(int f) {
    DiffractionImage img;
    img.metadata.width = kW;
    img.metadata.height = kH;
    img.pixels.resize(static_cast<std::size_t>(kW * kH));
    for (int r = 0; r < kH; ++r) {
        for (int c = 0; c < kW; ++c) {
            const bool shadow = c >= 140 && c < 180 && r >= 90 && r < 110;
            img.pixels[at(c, r)] = shadow ? 2 : 19 + (c + r + f) % 3;
        }
    }
    for (int c = 0; c < kW; ++c) {
        img.pixels[at(c, 0)] = kPixelMasked;
    }
    img.pixels[at(50, 60)] = kPixelOverloaded;                         // stuck
    img.pixels[at(70, 80)] = f % 2 == 0 ? kPixelOverloaded : 3;        // flickers
    img.pixels[at(72, 80)] = f % 4 == 1 ? kPixelOverloaded : 2;        // now and then
    img.pixels[at(90, 80)] = f == 3 || f == 4 ? kPixelOverloaded : f == 2 || f == 5 ? 4000 : 20;  // a spot
    img.pixels[at(120, 40)] = 5000;                                    // hot
    img.pixels[at(30, 30)] = f == 4 ? kPixelOverloaded : 20;           // one real spot
    img.pixels[at(31, 30)] = f == 4 ? 3000 : 20;
    return img;
}

mosflm::DetectorGeometry geometry() {
    mosflm::DetectorGeometry g;
    g.pixel_size_x_mm = 0.172;
    g.pixel_size_y_mm = 0.172;
    g.distance_mm = 100.0;
    g.wavelength_A = 1.0;
    g.nfast = kW;
    g.nslow = kH;
    return g;
}

/// d (Å) at pixel radius ``r_px`` for :func:`geometry`.
double d_at(double r_px) {
    return 1.0 / (2.0 * std::sin(0.5 * std::atan(r_px * 0.172 / 100.0)));
}

}  // namespace

int main() {
    using namespace mosflm::diffraction;

    // ----- regions (XDS pixels: 1-based, inclusive)
    {
        const auto rect = parse_untrusted_region({"RECTANGLE", "488", "494", "1", "2527"});
        CHECK(rect.contains(487, 0));    // XDS x 488
        CHECK(rect.contains(493, 2526)); // XDS x 494, y 2527
        CHECK(!rect.contains(486, 5));
        CHECK(!rect.contains(494, 5));
        const auto ell = parse_untrusted_region({"ellipse", "11", "31", "11", "21"});
        CHECK(ell.shape == UntrustedRegion::Shape::Ellipse);
        CHECK(ell.contains(20, 15));   // centre (21, 16)
        CHECK(!ell.contains(10, 10));  // bounding-box corner
        const auto quad = parse_untrusted_region({"QUAD", "1", "1", "11", "1", "11", "11", "1", "11"});
        CHECK(quad.contains(5, 5) && quad.contains(0, 0) && !quad.contains(11, 5));
        bool threw = false;
        try {
            (void)parse_untrusted_region({"RECTANGLE", "1", "2", "3"});
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        CHECK(threw);
        const ResolutionRange ring{3.93, 3.87};
        CHECK(ring.contains(3.9) && !ring.contains(3.8));
        CHECK_EQ(static_cast<int>(ice_ring_ranges().size()), 10);
    }

    // ----- sampling
    {
        const auto f = defpix_frames(1, 900, 10);
        CHECK_EQ(static_cast<int>(f.size()), 10);
        CHECK_EQ(f.front(), 1);
        CHECK_EQ(f.back(), 900);
        CHECK_EQ(static_cast<int>(defpix_frames(3, 5, 10).size()), 3);
    }

    // ----- a synthetic sweep, beam at the centre
    DefpixSettings settings;
    settings.untrusted.push_back(parse_untrusted_region({"RECTANGLE", "11", "20", "11", "20"}));
    settings.excluded_resolution.push_back({d_at(79.5), d_at(80.5)});
    std::vector<int> frames;
    for (int f = 0; f < 10; ++f) {
        frames.push_back(f);
    }
    const auto res = run_defpix(frame, frames, geometry(), 100 * 0.172, 100 * 0.172, settings);
    const auto status = [&](int c, int r) { return static_cast<PixelStatus>(res.status[at(c, r)]); };
    CHECK(status(5, 0) == PixelStatus::DetectorMasked);
    CHECK(status(50, 60) == PixelStatus::Defective);
    CHECK(status(70, 80) == PixelStatus::Defective);
    CHECK(status(72, 80) == PixelStatus::Defective);  // 3 of 10, background otherwise
    CHECK(status(90, 80) == PixelStatus::Trusted);    // 2 of 10, bright either side
    CHECK(status(120, 40) == PixelStatus::Hot);
    CHECK(status(30, 30) == PixelStatus::Trusted);  // saturated in one frame only
    CHECK(status(31, 30) == PixelStatus::Trusted);  // a spot in one frame only
    CHECK(status(15, 15) == PixelStatus::UntrustedRegion);
    CHECK_EQ(static_cast<int>(res.count(PixelStatus::UntrustedRegion)), 100);
    CHECK(status(180, 100) == PixelStatus::ExcludedResolution);
    CHECK(status(160, 100) == PixelStatus::Shaded);
    CHECK(status(100, 150) == PixelStatus::Trusted);
    CHECK_EQ(static_cast<int>(res.count(PixelStatus::Defective)), 3);
    CHECK_EQ(static_cast<int>(res.count(PixelStatus::Hot)), 1);
    CHECK(res.count(PixelStatus::Shaded) >= 40 * 16);
    CHECK(res.count(PixelStatus::Shaded) <= 48 * 24);
    CHECK_EQ(res.shadow_block_px, 8);
    // Background without each pixel's largest value: ~20 counts per frame.
    CHECK(std::fabs(res.background[at(100, 150)] - 20.0f) < 1.0f);
    const auto mask = res.trusted_mask();
    CHECK(mask[at(50, 60)] == 0 && mask[at(100, 150)] == 1 && mask[at(5, 0)] == 1);
    CHECK_EQ(res.n_rejected(), res.status.size() - res.count(PixelStatus::Trusted) -
                                   res.count(PixelStatus::DetectorMasked));

    // Two frames: only regions and resolution apply.
    {
        const auto two = run_defpix(frame, {0, 1}, geometry(), 100 * 0.172, 100 * 0.172, settings);
        CHECK_EQ(two.shadow_block_px, 0);
        CHECK(static_cast<PixelStatus>(two.status[at(50, 60)]) == PixelStatus::DetectorMasked);
        CHECK(static_cast<PixelStatus>(two.status[at(120, 40)]) == PixelStatus::Trusted);
        CHECK_EQ(static_cast<int>(two.count(PixelStatus::UntrustedRegion)), 100);
    }

    // ----- control image, and the CBF reader's pixel markers
    const auto dir = mosflm_test::output_dir("defpix");
    {
        write_defpix_cbf(dir / "defpix.cbf", res);
        const auto back = mosflm::image::CbfImageReader{}.read(dir / "defpix.cbf");
        CHECK_EQ(back.metadata.width, kW);
        CHECK_EQ(back.metadata.height, kH);
        CHECK(back.pixels[at(50, 60)] < 0 && back.pixels[at(5, 0)] < 0);
        CHECK(std::abs(back.pixels[at(100, 150)] - 200) <= 10);  // ~20 × 10 frames

        DiffractionImage img;
        img.metadata.width = 3;
        img.metadata.height = 2;
        img.pixels = {-2, -1, 999, 1000, 1048575, 5};
        mosflm::image::write_minicbf(dir / "cutoff.cbf", img,
                                     "# Pixel_size 172e-6 m x 172e-6 m\n# Count_cutoff 1000 counts\n");
        const auto read = mosflm::image::CbfImageReader{}.read(dir / "cutoff.cbf");
        CHECK_CLOSE(read.metadata.saturation_counts, 1000.0, 1e-9);
        const std::vector<std::int32_t> want{kPixelMasked, kPixelMasked, 999, kPixelOverloaded, kPixelOverloaded, 5};
        CHECK(read.pixels == want);
    }

    // ----- CONTROL keywords
    {
        mosflm::control::ControlSession s;
        s.parse_line("UNTRUSTED RECTANGLE 1 10 1 20");
        s.parse_line("UNTRUSTED ELLIPSE 100 140 100 120");
        s.parse_line("EXCLUDE RESOLUTION 3.93 3.87");
        s.parse_line("EXCLUDE ICE");
        s.parse_line("DEFPIX FRAMES 5 VALUE_RANGE 6000 25000 SHADOWS OFF FILE mask.cbf");
        const auto& d = s.state().defpix;
        CHECK_EQ(static_cast<int>(d.untrusted.size()), 2);
        CHECK_EQ(static_cast<int>(d.excluded_resolution.size()), 11);
        CHECK(d.enabled && d.frames == 5 && !d.find_shadows && d.find_defects);
        CHECK_CLOSE(d.value_range_low, 6000.0, 1e-9);
        CHECK(s.state().defpix_file == "mask.cbf");
        s.parse_line("DEFPIX OFF FILE OFF");
        s.parse_line("UNTRUSTED CLEAR");
        CHECK(!s.state().defpix.enabled && s.state().defpix_file.empty() && s.state().defpix.untrusted.empty());
        bool threw = false;
        try {
            s.parse_line("DEFPIX SIDEWAYS");
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        CHECK(threw);
    }

    // ----- XDS.INP import
    {
        const auto g = mosflm::io::XdsInpImporter::parse_string(
            "DATA_RANGE= 1 10\n"
            "UNTRUSTED_RECTANGLE= 488 494 1 2527\n"
            "UNTRUSTED_ELLIPSE= 1200 1290 1180 1250\n"
            "EXCLUDE_RESOLUTION_RANGE= 3.93 3.87\n"
            "VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS= 7000. 30000.\n"
        );
        const std::vector<std::string> want{"UNTRUSTED RECTANGLE 488 494 1 2527",
                                            "UNTRUSTED ELLIPSE 1200 1290 1180 1250",
                                            "EXCLUDE RESOLUTION 3.93 3.87",
                                            "DEFPIX VALUE_RANGE 7000. 30000."};
        CHECK(g.trusted_region_control == want);
        // XDS's "unknown" cell and space group are not a cell or space group.
        const auto unknown = mosflm::io::XdsInpImporter::parse_string(
            "DATA_RANGE= 1 10\nUNIT_CELL_CONSTANTS= 0 0 0 0 0 0\nSPACE_GROUP_NUMBER= 0\n");
        CHECK(!unknown.has_cell && !unknown.has_symmetry);
        mosflm::control::ControlSession s;
        for (const auto& line : g.to_control_lines()) {
            s.parse_line(line);
        }
        CHECK_EQ(static_cast<int>(s.state().defpix.untrusted.size()), 2);
        CHECK_CLOSE(s.state().defpix.value_range_high, 30000.0, 1e-9);
    }

    // ----- the beam halo: a flare filling the inner rings, so every ring is
    // its own median and the shaded test cannot see it
    {
        constexpr int kHw = 600;
        constexpr int kHc = 300;
        const auto halo_frame = [&](bool flare) {
            return [=](int f) {
                DiffractionImage img;
                img.metadata.width = kHw;
                img.metadata.height = kHw;
                img.pixels.resize(static_cast<std::size_t>(kHw * kHw));
                for (int r = 0; r < kHw; ++r) {
                    for (int c = 0; c < kHw; ++c) {
                        const double rad = std::hypot(c - kHc, r - kHc);
                        double v = 20.0 + (c + r + f) % 3;
                        if (rad <= 12.0) {
                            v = 1.0;  // the stop
                        } else if (flare) {
                            v += 60.0 * std::exp(-(rad - 12.0) / 25.0);
                        }
                        img.pixels[static_cast<std::size_t>(r) * kHw + static_cast<std::size_t>(c)] =
                            static_cast<std::int32_t>(std::lround(v));
                    }
                }
                return img;
            };
        };
        mosflm::DetectorGeometry g = geometry();
        g.nfast = kHw;
        g.nslow = kHw;
        DefpixSettings hs;
        hs.find_defects = false;
        const std::vector<int> hf{0, 1, 2, 3, 4};
        const auto with = run_defpix(halo_frame(true), hf, g, kHc * 0.172, kHc * 0.172, hs);
        CHECK(with.halo.found);
        // Excess 60 e^{-(r-12)/25} falls to 10 % of its peak near r ≈ 70 px.
        CHECK_MSG(with.halo.radius_px > 55.0 && with.halo.radius_px < 95.0,
                  "halo radius " + mosflm_test::fmt(with.halo.radius_px));
        const auto hstatus = [&](int c, int r) {
            return static_cast<PixelStatus>(with.status[static_cast<std::size_t>(r) * kHw + static_cast<std::size_t>(c)]);
        };
        CHECK(hstatus(kHc + 30, kHc) == PixelStatus::Halo);
        CHECK(hstatus(kHc, kHc + 150) == PixelStatus::Trusted);
        CHECK_EQ(static_cast<int>(with.trusted_mask()[static_cast<std::size_t>(kHc) * kHw + kHc + 30]), 0);
        // Spot finding may still look in the halo (against its local background).
        CHECK_EQ(static_cast<int>(with.spot_search_mask()[static_cast<std::size_t>(kHc) * kHw + kHc + 30]), 1);
        CHECK_CLOSE(with.halo_dmax_A, d_at(with.halo.radius_px + hs.halo_margin_px), 1e-6);

        const auto without = run_defpix(halo_frame(false), hf, g, kHc * 0.172, kHc * 0.172, hs);
        CHECK(!without.halo.found);
        CHECK_EQ(static_cast<int>(without.count(PixelStatus::Halo)), 0);
        CHECK(without.halo_dmax_A == 0.0);

        hs.find_halo = false;
        CHECK(!run_defpix(halo_frame(true), hf, g, kHc * 0.172, kHc * 0.172, hs).halo.found);
    }

    return mosflm_test::finish("defpix");
}
