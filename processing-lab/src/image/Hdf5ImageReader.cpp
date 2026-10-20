#include "mosflm/image/DiffractionImage.hpp"

#if !defined(MOSFLM_HAS_HDF5)

#include <stdexcept>

namespace mosflm::image {

DiffractionImage Hdf5ImageReader::read(
    const std::filesystem::path& path,
    int /*frame*/
) const {
    throw std::runtime_error(
        "ImageIO: HDF5/NeXus reader not built (need HDF5); cannot read " +
        path.string()
    );
}

int Hdf5ImageReader::num_images(const std::filesystem::path& path) const {
    throw std::runtime_error(
        "ImageIO: HDF5/NeXus reader not built (need HDF5); cannot read " +
        path.string()
    );
}

}  // namespace mosflm::image

#else

#include <hdf5.h>
#include <hdf5_hl.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

// Vendored filters (ProcessingLab bitshuffle/ + lz4/).
extern "C" {
#include "bshuf_h5filter.h"
extern const H5Z_class2_t H5Z_LZ4[1];
}

namespace mosflm::image {
namespace {

void ensure_filters_registered() {
    static std::once_flag once;
    std::call_once(once, [] {
        if (bshuf_register_h5filter() < 0) {
            // Fall through — H5Zregister may still succeed below.
        }
        if (H5Zregister(H5Z_LZ4) < 0) {
            // Non-fatal if already registered.
            H5Eclear2(H5E_DEFAULT);
        }
    });
}

/// Homebrew libhdf5 is built without thread safety (concurrent reads from
/// OpenMP FINDSPOTS / PROCESS crash), so every entry point holds this lock.
std::mutex& hdf5_mutex() {
    static std::mutex m;
    return m;
}

bool is_supported_hdf5_name(const std::filesystem::path& path) {
    const std::string name = path.filename().string();
    const std::string lower = [&] {
        std::string s = name;
        for (char& c : s) {
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        }
        return s;
    }();
    return lower.find("_master.h5") != std::string::npos ||
           (lower.size() >= 4 &&
            lower.compare(lower.size() - 4, 4, ".nxs") == 0) ||
           (lower.size() >= 5 &&
            lower.compare(lower.size() - 5, 5, ".hdf5") == 0) ||
           (lower.size() >= 3 &&
            lower.compare(lower.size() - 3, 3, ".h5") == 0);
}

hid_t open_file(const std::filesystem::path& path) {
    ensure_filters_registered();
    // Suppress HDF5 console noise for optional datasets.
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    const hid_t hdf = H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (hdf < 0) {
        throw std::runtime_error(
            "Hdf5ImageReader: failed to open " + path.string()
        );
    }
    return hdf;
}

bool read_int(hid_t loc, const char* path, int& out) {
    return H5LTread_dataset_int(loc, path, &out) >= 0;
}

bool read_double(hid_t loc, const char* path, double& out) {
    return H5LTread_dataset_double(loc, path, &out) >= 0;
}

bool read_string(hid_t loc, const char* path, std::string& out) {
    char buf[1024];
    std::memset(buf, 0, sizeof(buf));
    if (H5LTread_dataset_string(loc, path, buf) < 0) {
        return false;
    }
    out.assign(buf);
    // Trim trailing NULs already handled; strip trailing spaces.
    while (!out.empty() && (out.back() == '\0' || out.back() == ' ')) {
        out.pop_back();
    }
    return true;
}

/// NeXus ``vector`` attribute (3 doubles) on the object at ``obj_path``.
bool read_vector_attr(hid_t loc, const char* obj_path, double out[3]) {
    if (H5Aexists_by_name(loc, obj_path, "vector", H5P_DEFAULT) <= 0) {
        return false;
    }
    // get_attribute_info writes one dim per rank — check rank first.
    int rank = -1;
    if (H5LTget_attribute_ndims(loc, obj_path, "vector", &rank) < 0 ||
        rank != 1) {
        return false;
    }
    hsize_t dims[1] = {0};
    H5T_class_t cls{};
    std::size_t type_size = 0;
    if (H5LTget_attribute_info(
            loc, obj_path, "vector", dims, &cls, &type_size
        ) < 0 ||
        cls != H5T_FLOAT || dims[0] != 3) {
        return false;
    }
    return H5LTget_attribute_double(loc, obj_path, "vector", out) >= 0;
}

struct MasterMeta {
    int width{0};
    int height{0};
    double pixel_x_mm{0.075};
    double pixel_y_mm{0.075};
    double distance_mm{100.0};
    double wavelength_A{1.0};
    double beam_x_mm{0.0};
    double beam_y_mm{0.0};
    double osc_start_deg{0.0};
    double osc_range_deg{0.0};
    double rotation_axis_omega_deg{0.0};
    std::string detector_id{"Eiger"};
    int nimages{0};
    int bit_depth{16};
    unsigned int error_val{65535};
};

MasterMeta read_metadata(hid_t hdf, int frame_1based) {
    MasterMeta m;

    // --- detector size (classic Eiger vs DLS NXmx) -------------------------
    if (!read_int(
            hdf,
            "/entry/instrument/detector/detectorSpecific/x_pixels_in_detector",
            m.width
        )) {
        int size[2] = {0, 0};
        if (H5LTread_dataset_int(
                hdf, "/entry/instrument/detector/module/data_size", size
            ) >= 0) {
            // NXmx data_size is [fast, slow] in some files; thaumatin is
            // [4148, 4362] matching (fast, slow) with data dims (n, slow, fast).
            m.width = size[0];
            m.height = size[1];
        }
    }
    if (!read_int(
            hdf,
            "/entry/instrument/detector/detectorSpecific/y_pixels_in_detector",
            m.height
        ) &&
        m.height <= 0) {
        // already filled from data_size when present
    }

    // --- pixel size: classic stores metres; DLS NXmx often mm ---------------
    double px = -1.0;
    double py = -1.0;
    read_double(hdf, "/entry/instrument/detector/x_pixel_size", px);
    read_double(hdf, "/entry/instrument/detector/y_pixel_size", py);
    if (px > 0.0) {
        m.pixel_x_mm = (px < 0.01) ? px * 1000.0 : px;  // m→mm if SI
    }
    if (py > 0.0) {
        m.pixel_y_mm = (py < 0.01) ? py * 1000.0 : py;
    }

    // --- distance -----------------------------------------------------------
    double dist = -1.0;
    if (!read_double(hdf, "/entry/instrument/detector/distance", dist)) {
        read_double(hdf, "/entry/instrument/detector/detector_distance", dist);
    }
    if (dist > 0.0) {
        // Classic firmware: metres. If value > 10, treat as already mm (fromhdf5).
        m.distance_mm = (dist > 10.0) ? dist : dist * 1000.0;
    }

    // --- wavelength ---------------------------------------------------------
    if (!read_double(
            hdf, "/entry/sample/beam/incident_wavelength", m.wavelength_A
        )) {
        read_double(
            hdf, "/entry/instrument/beam/incident_wavelength", m.wavelength_A
        );
    }

    // --- beam centre (px) → mm ----------------------------------------------
    double beam_x_px = -1.0;
    double beam_y_px = -1.0;
    if (!read_double(hdf, "/entry/instrument/detector/beam_center_x", beam_x_px)) {
        read_double(hdf, "/entry/instrument/detector/beam_centre_x", beam_x_px);
    }
    if (!read_double(hdf, "/entry/instrument/detector/beam_center_y", beam_y_px)) {
        read_double(hdf, "/entry/instrument/detector/beam_centre_y", beam_y_px);
    }
    if (beam_x_px > 0.0 && m.pixel_x_mm > 0.0) {
        m.beam_x_mm = beam_x_px * m.pixel_x_mm;
    }
    if (beam_y_px > 0.0 && m.pixel_y_mm > 0.0) {
        m.beam_y_mm = beam_y_px * m.pixel_y_mm;
    }
    // Many NXmx masters omit beam_center_*; fall back to panel centre so
    // XDET/YDET and integration boxes are not pinned near (0,0).
    if (!(m.beam_x_mm > 0.0) && m.width > 0 && m.pixel_x_mm > 0.0) {
        m.beam_x_mm =
            0.5 * static_cast<double>(m.width) * m.pixel_x_mm;
    }
    if (!(m.beam_y_mm > 0.0) && m.height > 0 && m.pixel_y_mm > 0.0) {
        m.beam_y_mm =
            0.5 * static_cast<double>(m.height) * m.pixel_y_mm;
    }

    // --- bit depth / saturation ---------------------------------------------
    if (!read_int(hdf, "/entry/instrument/detector/bit_depth_image", m.bit_depth)) {
        m.bit_depth = 16;
    }
    if (m.bit_depth > 0 && m.bit_depth < 32) {
        m.error_val = static_cast<unsigned int>((1ull << m.bit_depth) - 1ull);
    }
    int sat = -1;
    if (read_int(hdf, "/entry/instrument/detector/saturation_value", sat) &&
        sat > 0) {
        m.error_val = static_cast<unsigned int>(sat);
    }

    std::string desc;
    if (read_string(hdf, "/entry/instrument/detector/description", desc)) {
        m.detector_id = desc;
    }

    // --- nimages ------------------------------------------------------------
    int nimages = 0;
    int ntrigger = 1;
    if (read_int(
            hdf,
            "/entry/instrument/detector/detectorSpecific/nimages",
            nimages
        )) {
        read_int(
            hdf,
            "/entry/instrument/detector/detectorSpecific/ntrigger",
            ntrigger
        );
        if (ntrigger > 1) {
            nimages *= ntrigger;
        }
        m.nimages = nimages;
    }

    // --- oscillation from omega/phi arrays ----------------------------------
    auto try_angle_array = [&](const char* path) -> bool {
        hid_t dset = H5Dopen2(hdf, path, H5P_DEFAULT);
        if (dset < 0) {
            return false;
        }
        hid_t space = H5Dget_space(dset);
        const int nd = H5Sget_simple_extent_ndims(space);
        hsize_t dims[1] = {0};
        if (nd != 1) {
            H5Sclose(space);
            H5Dclose(dset);
            return false;
        }
        H5Sget_simple_extent_dims(space, dims, nullptr);
        std::vector<double> angles(static_cast<std::size_t>(dims[0]));
        if (H5Dread(
                dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                angles.data()
            ) < 0) {
            H5Sclose(space);
            H5Dclose(dset);
            return false;
        }
        H5Sclose(space);
        H5Dclose(dset);
        if (angles.empty()) {
            return false;
        }
        const int idx = std::max(0, frame_1based - 1);
        if (idx >= static_cast<int>(angles.size())) {
            return false;
        }
        m.osc_start_deg = angles[static_cast<std::size_t>(idx)];
        if (idx + 1 < static_cast<int>(angles.size())) {
            m.osc_range_deg =
                angles[static_cast<std::size_t>(idx + 1)] - m.osc_start_deg;
        } else if (angles.size() >= 2) {
            m.osc_range_deg = angles[1] - angles[0];
        }
        if (m.nimages <= 0) {
            m.nimages = static_cast<int>(angles.size());
        }
        return true;
    };

    static const char* const kAnglePaths[] = {
        "/entry/data/omega",
        "/entry/sample/goniometer/omega",
        "/entry/sample/goniometer/phi",
        "/entry/sample/transformations/omega",
    };
    const char* angle_path = nullptr;
    for (const char* p : kAnglePaths) {
        if (try_angle_array(p)) {
            angle_path = p;
            break;
        }
    }

    // --- rotation axis vs detector fast/slow → MOSFLM detector omega --------
    // NeXus ``vector`` attributes share one lab frame, so only dot products
    // matter: ys (the φ axis in SpotToRlp) = −sin ω·x + cos ω·y must equal
    // the axis, giving ω = atan2(−axis·fast, axis·slow).
    double axis[3];
    double fast[3];
    double slow[3];
    const bool have_axis =
        (angle_path != nullptr && read_vector_attr(hdf, angle_path, axis)) ||
        read_vector_attr(hdf, "/entry/sample/transformations/omega", axis);
    if (have_axis &&
        read_vector_attr(
            hdf, "/entry/instrument/detector/module/fast_pixel_direction", fast
        ) &&
        read_vector_attr(
            hdf, "/entry/instrument/detector/module/slow_pixel_direction", slow
        )) {
        auto dot = [](const double* u, const double* v) {
            return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
        };
        const double norm = std::sqrt(dot(axis, axis) * dot(fast, fast));
        const double norm_s = std::sqrt(dot(axis, axis) * dot(slow, slow));
        if (norm > 1e-12 && norm_s > 1e-12) {
            const double af = dot(axis, fast) / norm;
            const double as = dot(axis, slow) / norm_s;
            // Axis roughly normal to the detector: ω is undefined, keep 0.
            if (std::hypot(af, as) > 0.5) {
                constexpr double kRad2Deg = 180.0 / 3.14159265358979323846;
                double om = std::atan2(-af, as) * kRad2Deg;
                if (om < 0.0) {
                    om += 360.0;
                }
                m.rotation_axis_omega_deg = om;
            }
        }
    }

    double osc_w = -1.0;
    if (read_double(
            hdf, "/entry/sample/goniometer/omega_range_average", osc_w
        ) &&
        osc_w > 0.0 && m.osc_range_deg <= 0.0) {
        m.osc_range_deg = osc_w;
    }

    // Infer size from first data block if still unknown.
    if (m.width <= 0 || m.height <= 0) {
        hid_t entry = H5Gopen2(hdf, "/entry", H5P_DEFAULT);
        hid_t group = H5Gopen2(entry, "data", H5P_DEFAULT);
        if (group < 0) {
            group = entry;
        }
        int block_start = 1;
        if (H5LTfind_dataset(group, "data_000000")) {
            block_start = 0;
        }
        char name[32];
        std::snprintf(name, sizeof(name), "data_%06d", block_start);
        hid_t data = H5Dopen2(group, name, H5P_DEFAULT);
        if (data >= 0) {
            hid_t space = H5Dget_space(data);
            hsize_t dims[3] = {0, 0, 0};
            if (H5Sget_simple_extent_ndims(space) == 3) {
                H5Sget_simple_extent_dims(space, dims, nullptr);
                m.height = static_cast<int>(dims[1]);
                m.width = static_cast<int>(dims[2]);
            }
            H5Sclose(space);
            H5Dclose(data);
        }
        if (group != entry) {
            H5Gclose(group);
        }
        H5Gclose(entry);
    }

    if (m.width <= 0 || m.height <= 0) {
        throw std::runtime_error(
            "Hdf5ImageReader: cannot determine detector size"
        );
    }
    return m;
}

struct DataLayout {
    hid_t group{-1};
    int block_start{1};
    int number_per_block{0};
    int total_frames{0};
};

DataLayout open_data_layout(hid_t hdf) {
    DataLayout layout;
    hid_t entry = H5Gopen2(hdf, "/entry", H5P_DEFAULT);
    if (entry < 0) {
        throw std::runtime_error("Hdf5ImageReader: /entry missing");
    }
    layout.group = H5Gopen2(entry, "data", H5P_DEFAULT);
    if (layout.group < 0) {
        layout.group = entry;
    } else {
        H5Gclose(entry);
    }

    if (H5LTfind_dataset(layout.group, "data_000000")) {
        layout.block_start = 0;
    } else {
        layout.block_start = 1;
    }

    char name[32];
    std::snprintf(name, sizeof(name), "data_%06d", layout.block_start);
    hid_t data = H5Dopen2(layout.group, name, H5P_DEFAULT);
    if (data < 0) {
        throw std::runtime_error(
            std::string("Hdf5ImageReader: failed to open /entry/data/") + name +
            " (are data_*.h5 files beside the master?)"
        );
    }
    hid_t space = H5Dget_space(data);
    if (H5Sget_simple_extent_ndims(space) != 3) {
        H5Sclose(space);
        H5Dclose(data);
        throw std::runtime_error(
            "Hdf5ImageReader: image dataset is not 3-D (n, slow, fast)"
        );
    }
    hsize_t dims[3] = {0, 0, 0};
    H5Sget_simple_extent_dims(space, dims, nullptr);
    layout.number_per_block = static_cast<int>(dims[0]);
    H5Sclose(space);
    H5Dclose(data);

    // Sum frames across data_* blocks.
    int block = layout.block_start;
    layout.total_frames = 0;
    for (;;) {
        std::snprintf(name, sizeof(name), "data_%06d", block);
        if (!H5LTfind_dataset(layout.group, name)) {
            break;
        }
        data = H5Dopen2(layout.group, name, H5P_DEFAULT);
        if (data < 0) {
            break;
        }
        space = H5Dget_space(data);
        H5Sget_simple_extent_dims(space, dims, nullptr);
        layout.total_frames += static_cast<int>(dims[0]);
        H5Sclose(space);
        H5Dclose(data);
        ++block;
    }
    return layout;
}

void read_frame_pixels(
    const DataLayout& layout,
    int frame_1based,
    int width,
    int height,
    unsigned int error_val,
    std::vector<std::int32_t>& pixels
) {
    if (frame_1based < 1) {
        throw std::runtime_error("Hdf5ImageReader: frame must be >= 1");
    }
    if (layout.number_per_block <= 0) {
        throw std::runtime_error("Hdf5ImageReader: empty data block");
    }
    if (layout.total_frames > 0 && frame_1based > layout.total_frames) {
        throw std::runtime_error(
            "Hdf5ImageReader: frame " + std::to_string(frame_1based) +
            " out of range (nimages=" + std::to_string(layout.total_frames) +
            ")"
        );
    }

    const int block_number =
        layout.block_start + (frame_1based - 1) / layout.number_per_block;
    const int frame_in_block =
        (frame_1based - 1) % layout.number_per_block;

    char name[32];
    std::snprintf(name, sizeof(name), "data_%06d", block_number);
    hid_t data = H5Dopen2(layout.group, name, H5P_DEFAULT);
    if (data < 0) {
        throw std::runtime_error(
            std::string("Hdf5ImageReader: failed to open ") + name
        );
    }
    hid_t dataspace = H5Dget_space(data);
    hsize_t dims[3] = {0, 0, 0};
    H5Sget_simple_extent_dims(dataspace, dims, nullptr);

    const hsize_t offset_in[3] = {
        static_cast<hsize_t>(frame_in_block), 0, 0
    };
    const hsize_t count[3] = {
        1, static_cast<hsize_t>(height), static_cast<hsize_t>(width)
    };
    const hsize_t mem_dims[3] = {
        1, static_cast<hsize_t>(height), static_cast<hsize_t>(width)
    };
    hid_t memspace = H5Screate_simple(3, mem_dims, nullptr);
    if (H5Sselect_hyperslab(
            dataspace, H5S_SELECT_SET, offset_in, nullptr, count, nullptr
        ) < 0) {
        H5Sclose(memspace);
        H5Sclose(dataspace);
        H5Dclose(data);
        throw std::runtime_error("Hdf5ImageReader: hyperslab select failed");
    }
    const hsize_t offset_out[3] = {0, 0, 0};
    H5Sselect_hyperslab(
        memspace, H5S_SELECT_SET, offset_out, nullptr, count, nullptr
    );

    const std::size_t npix =
        static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    std::vector<unsigned int> buf(npix);
    if (H5Dread(
            data, H5T_NATIVE_UINT, memspace, dataspace, H5P_DEFAULT, buf.data()
        ) < 0) {
        H5Sclose(memspace);
        H5Sclose(dataspace);
        H5Dclose(data);
        throw std::runtime_error(
            "Hdf5ImageReader: H5Dread failed (filters / frame / data files?)"
        );
    }
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(data);

    pixels.resize(npix);
    for (std::size_t i = 0; i < npix; ++i) {
        if (buf[i] == error_val || buf[i] == static_cast<unsigned int>(-1)) {
            pixels[i] = -1;
        } else {
            pixels[i] = static_cast<std::int32_t>(buf[i]);
        }
    }
}

}  // namespace

DiffractionImage Hdf5ImageReader::read(
    const std::filesystem::path& path,
    int frame
) const {
    if (!is_supported_hdf5_name(path)) {
        throw std::runtime_error(
            "Hdf5ImageReader: expected *_master.h5 or .nxs, got " + path.string()
        );
    }
    const std::lock_guard<std::mutex> lock(hdf5_mutex());
    const hid_t hdf = open_file(path);
    try {
        MasterMeta meta = read_metadata(hdf, frame);
        DataLayout layout = open_data_layout(hdf);
        if (meta.nimages <= 0) {
            meta.nimages = layout.total_frames;
        }

        DiffractionImage out;
        out.metadata.width = meta.width;
        out.metadata.height = meta.height;
        out.metadata.pixel_size_x_mm = meta.pixel_x_mm;
        out.metadata.pixel_size_y_mm = meta.pixel_y_mm;
        out.metadata.detector_distance_mm = meta.distance_mm;
        out.metadata.wavelength_angstrom = meta.wavelength_A;
        out.metadata.beam_x_mm = meta.beam_x_mm;
        out.metadata.beam_y_mm = meta.beam_y_mm;
        out.metadata.oscillation_start_deg = meta.osc_start_deg;
        out.metadata.oscillation_range_deg =
            meta.osc_range_deg > 0.0 ? meta.osc_range_deg : 0.0;
        out.metadata.rotation_axis_omega_deg = meta.rotation_axis_omega_deg;
        out.metadata.detector_id = meta.detector_id;

        read_frame_pixels(
            layout, frame, meta.width, meta.height, meta.error_val, out.pixels
        );

        H5Gclose(layout.group);
        H5Fclose(hdf);
        return out;
    } catch (...) {
        H5Fclose(hdf);
        throw;
    }
}

int Hdf5ImageReader::num_images(const std::filesystem::path& path) const {
    if (!is_supported_hdf5_name(path)) {
        throw std::runtime_error(
            "Hdf5ImageReader: expected *_master.h5 or .nxs, got " + path.string()
        );
    }
    const std::lock_guard<std::mutex> lock(hdf5_mutex());
    const hid_t hdf = open_file(path);
    try {
        DataLayout layout = open_data_layout(hdf);
        const int n = layout.total_frames;
        H5Gclose(layout.group);
        H5Fclose(hdf);
        return n;
    } catch (...) {
        H5Fclose(hdf);
        throw;
    }
}

}  // namespace mosflm::image

#endif  // MOSFLM_HAS_HDF5
