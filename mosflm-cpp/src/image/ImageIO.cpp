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

bool is_hdf5_ext(const std::string& ext) {
    return ext == ".h5" || ext == ".nxs" || ext == ".hdf5";
}

}  // namespace

DiffractionImage read_diffraction_image(
    const std::filesystem::path& path,
    int frame
) {
    const std::string ext = lower_ext(path);
    if (ext == ".cbf") {
        return CbfImageReader{}.read(path);
    }
    if (is_hdf5_ext(ext)) {
#if defined(MOSFLM_HAS_HDF5)
        return Hdf5ImageReader{}.read(path, frame);
#else
        throw std::runtime_error(
            "ImageIO: HDF5/NeXus reader not built for " + path.string() +
            " (install HDF5 and rebuild with MOSFLM_HAS_HDF5)"
        );
#endif
    }
    if (ext == ".jpg" || ext == ".jpeg") {
        throw std::runtime_error(
            "ImageIO: JPEG diffraction frames not supported (" + path.string() +
            ")"
        );
    }
    // Default: ADSC/SMV (``.img``, ``.smv``, or extension-less CCP4 examples).
    (void)frame;
    return SmvImageReader{}.read(path);
}

}  // namespace mosflm::image
