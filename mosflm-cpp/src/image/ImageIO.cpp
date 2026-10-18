#include "mosflm/image/ImageIO.hpp"

#include <cctype>
#include <stdexcept>
#include <string>

namespace mosflm::image {
namespace {

std::string lower_ext(const std::filesystem::path& path) {
    std::string e = path.extension().string();
    for (char& c : e) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return e;
}

}  // namespace

DiffractionImage read_diffraction_image(const std::filesystem::path& path) {
    const std::string ext = lower_ext(path);
    if (ext == ".cbf") {
        return CbfImageReader{}.read(path);
    }
    if (ext == ".h5" || ext == ".nxs" || ext == ".hdf5") {
        throw std::runtime_error(
            "ImageIO: HDF5/NeXus reader not yet implemented for " + path.string() +
            " (Phase E / E5)"
        );
    }
    if (ext == ".jpg" || ext == ".jpeg") {
        throw std::runtime_error(
            "ImageIO: JPEG diffraction frames not supported (" + path.string() +
            ")"
        );
    }
    // Default: ADSC/SMV (``.img``, ``.smv``, or extension-less CCP4 examples).
    return SmvImageReader{}.read(path);
}

}  // namespace mosflm::image
