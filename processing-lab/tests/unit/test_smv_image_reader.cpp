#include "mosflm/image/DiffractionImage.hpp"
#include "test_util.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>

using mosflm::image::DiffractionImage;
using mosflm::image::SmvImageReader;

namespace {

std::filesystem::path find_insulin_img() {
    if (const char* env = std::getenv("MOSFLM_IMAGE")) {
        return env;
    }
    if (const char* ccp4 = std::getenv("CCP4")) {
        const auto p =
            std::filesystem::path(ccp4) / "examples/data/insulin_1_001.img";
        if (std::filesystem::exists(p)) {
            return p;
        }
    }
    const std::filesystem::path fallback =
        "/Applications/ccp4-9/examples/data/insulin_1_001.img";
    return fallback;
}

}  // namespace

int main() {
    const auto path = find_insulin_img();
    CHECK_MSG(std::filesystem::exists(path), "missing insulin SMV (set MOSFLM_IMAGE or CCP4)");

    const DiffractionImage image = SmvImageReader{}.read(path);
    CHECK_EQ(image.metadata.width, 2304);
    CHECK_EQ(image.metadata.height, 2304);
    CHECK_EQ(
        static_cast<int>(image.pixels.size()), 2304 * 2304
    );
    CHECK_CLOSE(image.metadata.pixel_size_x_mm, 0.0816, 1e-6);
    CHECK_CLOSE(image.metadata.detector_distance_mm, 160.0, 1e-3);
    CHECK_CLOSE(image.metadata.wavelength_angstrom, 0.9790, 1e-4);
    CHECK_CLOSE(image.metadata.beam_x_mm, 94.0, 1e-3);
    CHECK_CLOSE(image.metadata.beam_y_mm, 94.7, 1e-3);
    CHECK_CLOSE(image.metadata.oscillation_start_deg, 0.0, 1e-6);
    CHECK_CLOSE(image.metadata.oscillation_range_deg, 1.0, 1e-6);

    // Beam-stop / edges may be zero; expect real counts elsewhere.
    std::int32_t max_px = 0;
    std::int64_t sum = 0;
    for (std::int32_t v : image.pixels) {
        max_px = std::max(max_px, v);
        sum += v;
    }
    CHECK(max_px > 1000);
    CHECK(max_px < 65535);
    const double mean =
        static_cast<double>(sum) / static_cast<double>(image.pixels.size());
    CHECK(mean > 50.0);

    const auto det = image.metadata.to_detector_geometry();
    CHECK_EQ(det.nfast, 2304);
    CHECK_EQ(det.nslow, 2304);

    std::cout << "smv_image_reader: " << path << " ok\n";
    return mosflm_test::finish("smv_image_reader");
}
