#include "mosflm/image/ImageIO.hpp"
#include "test_util.hpp"

#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>

using mosflm::image::read_diffraction_image;

int main() {
    bool threw = false;
    try {
        (void)read_diffraction_image("missing.h5");
    } catch (const std::runtime_error& ex) {
        threw = true;
        CHECK_MSG(std::string(ex.what()).find("HDF5") != std::string::npos,
                  "expected HDF5 message");
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

    return mosflm_test::finish("image_io");
}
