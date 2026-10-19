#include "mosflm/image/DiffractionImage.hpp"
#include "test_util.hpp"

#include <filesystem>
#include <string>

using mosflm::image::CbfImageReader;
using mosflm::image::DiffractionImage;

int main() {
    const char* env = std::getenv("MOSFLM_CBF_FIXTURE");
    if (env == nullptr) {
        std::cerr << "MOSFLM_CBF_FIXTURE not set\n";
        return 1;
    }
    const std::filesystem::path fixture(env);
    CHECK(std::filesystem::exists(fixture));

    const DiffractionImage image = CbfImageReader{}.read(fixture);

    CHECK_EQ(image.metadata.width, 128);
    CHECK_EQ(image.metadata.height, 128);
    CHECK_EQ(static_cast<int>(image.pixels.size()), 128 * 128);

    CHECK_CLOSE(image.metadata.wavelength_angstrom, 0.979, 1e-6);
    CHECK_CLOSE(image.metadata.pixel_size_x_mm, 0.172, 1e-6);
    CHECK_CLOSE(image.metadata.pixel_size_y_mm, 0.172, 1e-6);
    CHECK_CLOSE(image.metadata.detector_distance_mm, 160.0, 1e-3);
    CHECK_CLOSE(image.metadata.beam_x_mm, 64.0 * 0.172, 1e-3);
    CHECK_CLOSE(image.metadata.beam_y_mm, 64.0 * 0.172, 1e-3);
    CHECK_CLOSE(image.metadata.oscillation_start_deg, 0.0, 1e-9);
    CHECK_CLOSE(image.metadata.oscillation_range_deg, 1.0, 1e-9);
    CHECK(image.metadata.detector_id == "SYNTHETIC-MINI");

    // Background ~50, plus the planted peak at the centre.
    CHECK_EQ(static_cast<int>(image.pixels[0]), 50);
    const std::size_t centre =
        static_cast<std::size_t>(64) * 128u + static_cast<std::size_t>(64);
    CHECK(image.pixels[centre] > 500);

    // Detector geometry mapping roundtrip.
    const auto det = image.metadata.to_detector_geometry();
    CHECK_EQ(det.nfast, 128);
    CHECK_EQ(det.nslow, 128);
    CHECK_CLOSE(det.wavelength_A, 0.979, 1e-6);

    return mosflm_test::finish("cbf_image_reader");
}
