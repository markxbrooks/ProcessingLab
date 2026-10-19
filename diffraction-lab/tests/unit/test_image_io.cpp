#include "mosflm/image/ImageIO.hpp"
#include "test_util.hpp"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

using mosflm::image::read_diffraction_image;

namespace {

std::filesystem::path find_hdf5_fixture() {
    if (const char* env = std::getenv("MOSFLM_HDF5_FIXTURE")) {
        return env;
    }
    const char* candidates[] = {
        "/Users/brooks/Downloads/thaumatin_9_1_master.h5",
    };
    for (const char* c : candidates) {
        if (std::filesystem::exists(c)) {
            return c;
        }
    }
    return {};
}

}  // namespace

int main() {
    bool threw = false;
    try {
        (void)read_diffraction_image("missing.h5");
    } catch (const std::runtime_error& ex) {
        threw = true;
        CHECK_MSG(std::string(ex.what()).find("HDF5") != std::string::npos ||
                      std::string(ex.what()).find("Hdf5") != std::string::npos ||
                      std::string(ex.what()).find("failed to open") !=
                          std::string::npos,
                  "expected HDF5-related error");
    }
    CHECK(threw);

    const char* ccp4 = std::getenv("CCP4");
    const std::filesystem::path img =
        ccp4 ? std::filesystem::path(ccp4) / "examples/data/insulin_1_001.img"
             : std::filesystem::path(
                   "/Applications/ccp4-9/examples/data/insulin_1_001.img"
               );
    if (std::filesystem::exists(img)) {
        const auto image = read_diffraction_image(img);
        CHECK(image.metadata.width > 0);
        CHECK(image.metadata.height > 0);
    }

#if defined(MOSFLM_HAS_HDF5)
    const auto h5 = find_hdf5_fixture();
    if (h5.empty() || !std::filesystem::exists(h5)) {
        std::cout << "image_io: SKIP HDF5 fixture (set MOSFLM_HDF5_FIXTURE)\n";
    } else {
        const auto image = read_diffraction_image(h5, 1);
        CHECK_MSG(image.metadata.width > 100, "HDF5 width");
        CHECK_MSG(image.metadata.height > 100, "HDF5 height");
        CHECK_MSG(
            image.pixels.size() ==
                static_cast<std::size_t>(image.metadata.width) *
                    static_cast<std::size_t>(image.metadata.height),
            "HDF5 pixel count"
        );
        CHECK_MSG(image.metadata.pixel_size_x_mm > 0.0, "HDF5 pixel size");
        CHECK_MSG(
            image.metadata.detector_distance_mm > 0.0, "HDF5 distance"
        );
        CHECK_MSG(
            image.metadata.wavelength_angstrom > 0.0, "HDF5 wavelength"
        );
        CHECK_MSG(image.metadata.beam_x_mm > 0.0, "HDF5 beam x");
        CHECK_MSG(image.metadata.beam_y_mm > 0.0, "HDF5 beam y");
        std::cout << "image_io: HDF5 fixture=" << h5
                  << " size=" << image.metadata.width << "x"
                  << image.metadata.height
                  << " dist=" << image.metadata.detector_distance_mm
                  << " wl=" << image.metadata.wavelength_angstrom
                  << " beam=(" << image.metadata.beam_x_mm << ","
                  << image.metadata.beam_y_mm << ")"
                  << " osc=" << image.metadata.oscillation_start_deg << "+"
                  << image.metadata.oscillation_range_deg << '\n';
    }
#else
    std::cout << "image_io: SKIP HDF5 (not built)\n";
#endif

    return mosflm_test::finish("image_io");
}
