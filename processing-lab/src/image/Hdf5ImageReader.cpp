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

// Portions of this file (NXdetector / NXdata discovery, the single NXmx
// ``data`` dataset, ``units`` scaling, ``pixel_mask`` handling and the direct
// bitshuffle-LZ4 chunk read) are adapted from durin, the XDS HDF5 plugin:
//
//   Copyright (c) 2018-2019, Diamond Light Source Ltd. All rights reserved.
//   Author: Charles Mita
//
// durin is distributed under the BSD 3-Clause licence; its copyright notice,
// conditions and disclaimer are reproduced in third_party/durin/LICENSE in
// this repository.

#include <hdf5.h>
#include <hdf5_hl.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

// Vendored filters (third_party/bitshuffle/ + third_party/lz4/).
extern "C" {
#include "bitshuffle.h"
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
/// OpenMP FINDSPOTS / PROCESS crash), so every HDF5 call holds this lock.
/// Bitshuffle-LZ4 decoding and pixel conversion run after it is released.
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

bool read_int(hid_t loc, const std::string& path, int& out) {
    return H5LTread_dataset_int(loc, path.c_str(), &out) >= 0;
}

bool read_double(hid_t loc, const std::string& path, double& out) {
    return H5LTread_dataset_double(loc, path.c_str(), &out) >= 0;
}

bool read_string(hid_t loc, const std::string& path, std::string& out) {
    char buf[1024];
    std::memset(buf, 0, sizeof(buf));
    if (H5LTread_dataset_string(loc, path.c_str(), buf) < 0) {
        return false;
    }
    out.assign(buf);
    // Trim trailing NULs already handled; strip trailing spaces.
    while (!out.empty() && (out.back() == '\0' || out.back() == ' ')) {
        out.pop_back();
    }
    return true;
}

/// Scalar string attribute, fixed or variable length (durin read_pixel_info).
bool read_string_attr(hid_t obj, const char* name, std::string& out) {
    if (H5Aexists(obj, name) <= 0) {
        return false;
    }
    const hid_t attr = H5Aopen(obj, name, H5P_DEFAULT);
    if (attr < 0) {
        return false;
    }
    bool ok = false;
    const hid_t space = H5Aget_space(attr);
    const hid_t type = H5Aget_type(attr);
    if (space >= 0 && type >= 0 && H5Tget_class(type) == H5T_STRING &&
        H5Sget_simple_extent_npoints(space) == 1) {
        const hid_t mem = H5Tcopy(H5T_C_S1);
        if (H5Tis_variable_str(type) > 0) {
            H5Tset_size(mem, H5T_VARIABLE);
            char* s = nullptr;
            if (H5Aread(attr, mem, &s) >= 0 && s != nullptr) {
                out.assign(s);
                H5free_memory(s);
                ok = true;
            }
        } else {
            // Do not assume the file left room for a terminating NUL.
            const std::size_t n = H5Tget_size(type);
            std::vector<char> buf(n + 1, '\0');
            H5Tset_size(mem, n);
            if (H5Aread(attr, mem, buf.data()) >= 0) {
                out.assign(buf.data());
                ok = true;
            }
        }
        H5Tclose(mem);
    }
    if (type >= 0) {
        H5Tclose(type);
    }
    if (space >= 0) {
        H5Sclose(space);
    }
    H5Aclose(attr);
    while (!out.empty() && out.back() == ' ') {
        out.pop_back();
    }
    return ok;
}

/// Millimetres per ``units`` string, or 0 if unrecognised (durin
/// scale_from_units, which returns metres).
double mm_per_unit(std::string units) {
    for (char& c : units) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    if (units == "m" || units == "metre" || units == "metres" ||
        units == "meter" || units == "meters") {
        return 1000.0;
    }
    if (units == "cm" || units == "centimetres" || units == "centimeters") {
        return 10.0;
    }
    if (units == "mm" || units == "millimetres" || units == "millimeters") {
        return 1.0;
    }
    if (units == "um" || units == "microns" || units == "micrometres" ||
        units == "micrometers") {
        return 0.001;
    }
    return 0.0;
}

/// Length dataset in mm, scaled by its ``units`` attribute. Without a
/// recognised unit, ``guess_mm`` converts the bare value.
bool read_length_mm(
    hid_t loc,
    const std::string& path,
    double& out_mm,
    double (*guess_mm)(double)
) {
    double value = 0.0;
    if (!read_double(loc, path, value)) {
        return false;
    }
    double scale = 0.0;
    const hid_t dset = H5Dopen2(loc, path.c_str(), H5P_DEFAULT);
    if (dset >= 0) {
        std::string units;
        if (read_string_attr(dset, "units", units)) {
            scale = mm_per_unit(units);
        }
        H5Dclose(dset);
    }
    out_mm = scale > 0.0 ? value * scale : guess_mm(value);
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

/// True if ``path`` names a group (resolving soft / external links).
bool is_group(hid_t hdf, const std::string& path) {
    if (H5Lexists(hdf, path.c_str(), H5P_DEFAULT) <= 0) {
        return false;
    }
    const hid_t g = H5Gopen2(hdf, path.c_str(), H5P_DEFAULT);
    if (g < 0) {
        return false;
    }
    H5Gclose(g);
    return true;
}

/// ``preferred`` if it is a group, else the first child of ``parent`` whose
/// ``NX_class`` is ``nx_class`` (durin det_visit_callback). Empty if none.
std::string find_nx_group(
    hid_t hdf,
    const std::string& parent,
    const std::string& preferred,
    const char* nx_class
) {
    if (is_group(hdf, preferred)) {
        return preferred;
    }
    if (!is_group(hdf, parent)) {
        return {};
    }
    const hid_t g = H5Gopen2(hdf, parent.c_str(), H5P_DEFAULT);
    H5G_info_t info{};
    std::string found;
    if (H5Gget_info(g, &info) >= 0) {
        for (hsize_t i = 0; i < info.nlinks && found.empty(); ++i) {
            char name[256] = {0};
            if (H5Lget_name_by_idx(
                    g, ".", H5_INDEX_NAME, H5_ITER_INC, i, name, sizeof(name),
                    H5P_DEFAULT
                ) < 0) {
                continue;
            }
            const std::string child = parent + "/" + name;
            const hid_t c = H5Gopen2(hdf, child.c_str(), H5P_DEFAULT);
            if (c < 0) {
                continue;
            }
            std::string cls;
            if (read_string_attr(c, "NX_class", cls) && cls == nx_class) {
                found = child;
            }
            H5Gclose(c);
        }
    }
    H5Gclose(g);
    return found;
}

/// One 3-D (n, slow, fast) image dataset in the data group.
struct DataBlock {
    std::string name;
    int frames{0};
};

/// DECTRIS ``data_NNNNNN`` blocks, or the single NXmx ``data`` dataset
/// (often a virtual dataset over the data files).
struct DataLayout {
    hid_t group{-1};
    std::vector<DataBlock> blocks;
    int total_frames{0};
    int width{0};
    int height{0};
};

/// Frames in a 3-D dataset, filling ``height`` / ``width``; −1 if not 3-D or
/// unreadable (e.g. a missing external data file).
int dataset_frames(hid_t group, const char* name, int& height, int& width) {
    const hid_t data = H5Dopen2(group, name, H5P_DEFAULT);
    if (data < 0) {
        return -1;
    }
    const hid_t space = H5Dget_space(data);
    hsize_t dims[3] = {0, 0, 0};
    int frames = -1;
    if (H5Sget_simple_extent_ndims(space) == 3) {
        H5Sget_simple_extent_dims(space, dims, nullptr);
        frames = static_cast<int>(dims[0]);
        height = static_cast<int>(dims[1]);
        width = static_cast<int>(dims[2]);
    }
    H5Sclose(space);
    H5Dclose(data);
    return frames;
}

DataLayout open_data_layout(hid_t hdf) {
    DataLayout layout;
    if (!is_group(hdf, "/entry")) {
        throw std::runtime_error("Hdf5ImageReader: /entry missing");
    }
    std::string data_path = find_nx_group(hdf, "/entry", "/entry/data", "NXdata");
    if (data_path.empty()) {
        data_path = "/entry";
    }
    layout.group = H5Gopen2(hdf, data_path.c_str(), H5P_DEFAULT);

    const int block_start =
        H5LTfind_dataset(layout.group, "data_000000") ? 0 : 1;
    char name[32];
    std::snprintf(name, sizeof(name), "data_%06d", block_start);
    if (H5LTfind_dataset(layout.group, name)) {
        for (int block = block_start;; ++block) {
            std::snprintf(name, sizeof(name), "data_%06d", block);
            if (!H5LTfind_dataset(layout.group, name)) {
                break;
            }
            const int n =
                dataset_frames(layout.group, name, layout.height, layout.width);
            if (n < 0) {
                if (block == block_start) {
                    H5Gclose(layout.group);
                    throw std::runtime_error(
                        "Hdf5ImageReader: failed to open " + data_path + "/" +
                        name + " (are data_*.h5 files beside the master?)"
                    );
                }
                break;
            }
            layout.blocks.push_back({name, n});
            layout.total_frames += n;
        }
    } else if (H5LTfind_dataset(layout.group, "data")) {
        const int n =
            dataset_frames(layout.group, "data", layout.height, layout.width);
        if (n < 0) {
            H5Gclose(layout.group);
            throw std::runtime_error(
                "Hdf5ImageReader: " + data_path +
                "/data is not a readable 3-D (n, slow, fast) image dataset"
            );
        }
        layout.blocks.push_back({"data", n});
        layout.total_frames = n;
    } else {
        H5Gclose(layout.group);
        throw std::runtime_error(
            "Hdf5ImageReader: no data_NNNNNN or data image dataset in " +
            data_path
        );
    }
    return layout;
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

MasterMeta read_metadata(
    hid_t hdf,
    const std::string& det,
    int frame_1based,
    const DataLayout& layout
) {
    MasterMeta m;

    // --- detector size (classic Eiger vs DLS NXmx) -------------------------
    if (!read_int(hdf, det + "/detectorSpecific/x_pixels_in_detector", m.width)) {
        int size[2] = {0, 0};
        if (H5LTread_dataset_int(
                hdf, (det + "/module/data_size").c_str(), size
            ) >= 0) {
            // NXmx data_size is [fast, slow] in some files; thaumatin is
            // [4148, 4362] matching (fast, slow) with data dims (n, slow, fast).
            m.width = size[0];
            m.height = size[1];
        }
    }
    if (!read_int(hdf, det + "/detectorSpecific/y_pixels_in_detector", m.height) &&
        m.height <= 0) {
        // already filled from data_size when present
    }

    // --- pixel size: ``units`` when present; else classic stores metres and
    // DLS NXmx often mm -------------------------------------------------------
    double px = -1.0;
    double py = -1.0;
    const auto pixel_guess = [](double v) { return v < 0.01 ? v * 1000.0 : v; };
    read_length_mm(hdf, det + "/x_pixel_size", px, pixel_guess);
    read_length_mm(hdf, det + "/y_pixel_size", py, pixel_guess);
    if (px > 0.0) {
        m.pixel_x_mm = px;
    }
    if (py > 0.0) {
        m.pixel_y_mm = py;
    }

    // --- distance -----------------------------------------------------------
    // Without ``units``: classic firmware stores metres; a value > 10 is
    // taken as already mm (fromhdf5).
    double dist = -1.0;
    const auto dist_guess = [](double v) { return v > 10.0 ? v : v * 1000.0; };
    if (!read_length_mm(hdf, det + "/distance", dist, dist_guess)) {
        read_length_mm(hdf, det + "/detector_distance", dist, dist_guess);
    }
    if (dist > 0.0) {
        m.distance_mm = dist;
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
    if (!read_double(hdf, det + "/beam_center_x", beam_x_px)) {
        read_double(hdf, det + "/beam_centre_x", beam_x_px);
    }
    if (!read_double(hdf, det + "/beam_center_y", beam_y_px)) {
        read_double(hdf, det + "/beam_centre_y", beam_y_px);
    }
    if (beam_x_px > 0.0 && m.pixel_x_mm > 0.0) {
        m.beam_x_mm = beam_x_px * m.pixel_x_mm;
    }
    if (beam_y_px > 0.0 && m.pixel_y_mm > 0.0) {
        m.beam_y_mm = beam_y_px * m.pixel_y_mm;
    }

    // Infer size from the image datasets if still unknown.
    if (m.width <= 0 || m.height <= 0) {
        m.width = layout.width;
        m.height = layout.height;
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
    if (!read_int(hdf, det + "/bit_depth_image", m.bit_depth)) {
        m.bit_depth = 16;
    }
    if (m.bit_depth > 0 && m.bit_depth < 32) {
        m.error_val = static_cast<unsigned int>((1ull << m.bit_depth) - 1ull);
    }
    int sat = -1;
    if (read_int(hdf, det + "/saturation_value", sat) && sat > 0) {
        m.error_val = static_cast<unsigned int>(sat);
    }

    std::string desc;
    if (read_string(hdf, det + "/description", desc)) {
        m.detector_id = desc;
    }

    // --- nimages ------------------------------------------------------------
    int nimages = 0;
    int ntrigger = 1;
    if (read_int(hdf, det + "/detectorSpecific/nimages", nimages)) {
        read_int(hdf, det + "/detectorSpecific/ntrigger", ntrigger);
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
            hdf, (det + "/module/fast_pixel_direction").c_str(), fast
        ) &&
        read_vector_attr(
            hdf, (det + "/module/slow_pixel_direction").c_str(), slow
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

    if (m.width <= 0 || m.height <= 0) {
        throw std::runtime_error(
            "Hdf5ImageReader: cannot determine detector size"
        );
    }
    return m;
}

/// Per-pixel "ignore" flags from ``pixel_mask``; null when the file has none.
using PixelMask = std::vector<std::uint8_t>;

/// NXmx ``pixel_mask`` or DECTRIS ``detectorSpecific/pixel_mask`` (durin
/// get_nxs_pixel_mask / get_dectris_eiger_pixel_mask). Any bit in the low
/// byte (gap, dead, cold, hot, noisy, …) marks the pixel as not to be used.
/// Cached per file; call with ``hdf5_mutex`` held.
std::shared_ptr<const PixelMask> read_pixel_mask(
    hid_t hdf,
    const std::filesystem::path& path,
    const std::string& det,
    int width,
    int height
) {
    static std::map<std::string, std::shared_ptr<const PixelMask>> cache;
    std::error_code ec;
    const auto mtime = std::filesystem::last_write_time(path, ec);
    const std::string key =
        path.string() + '|' +
        std::to_string(ec ? 0LL : static_cast<long long>(
                                     mtime.time_since_epoch().count()));
    if (const auto it = cache.find(key); it != cache.end()) {
        return it->second;
    }

    std::shared_ptr<const PixelMask> result;
    for (const std::string p :
         {det + "/pixel_mask", det + "/detectorSpecific/pixel_mask"}) {
        // A dangling external link (mask file not written) exists but
        // does not open; try the next location.
        if (H5Lexists(hdf, p.c_str(), H5P_DEFAULT) <= 0) {
            continue;
        }
        const hid_t dset = H5Dopen2(hdf, p.c_str(), H5P_DEFAULT);
        if (dset < 0) {
            continue;
        }
        const hid_t space = H5Dget_space(dset);
        hsize_t dims[2] = {0, 0};
        const bool shape_ok =
            H5Sget_simple_extent_ndims(space) == 2 &&
            H5Sget_simple_extent_dims(space, dims, nullptr) == 2 &&
            dims[0] == static_cast<hsize_t>(height) &&
            dims[1] == static_cast<hsize_t>(width);
        H5Sclose(space);
        if (shape_ok) {
            const std::size_t npix =
                static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
            std::vector<std::uint32_t> bits(npix);
            if (H5Dread(
                    dset, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    bits.data()
                ) >= 0) {
                auto mask = std::make_shared<PixelMask>(npix);
                for (std::size_t i = 0; i < npix; ++i) {
                    (*mask)[i] = (bits[i] & UINT32_C(0x000000FF)) != 0;
                }
                result = std::move(mask);
            }
        }
        H5Dclose(dset);
        if (result) {
            break;
        }
    }
    cache.emplace(key, result);
    return result;
}

/// One frame as read from the file, before decoding. Either ``values``
/// (already converted by H5Dread) or ``chunk`` (raw chunk bytes) is filled.
struct RawFrame {
    std::vector<unsigned int> values;
    std::vector<unsigned char> chunk;
    bool lz4{false};          ///< chunk is bitshuffle-LZ4 compressed
    std::size_t elem_size{0};  ///< bytes per pixel in the file
    std::size_t npix{0};
};

std::uint64_t read_uint64_be(const unsigned char* p) {
    std::uint64_t v = 0;
    for (int i = 0; i < 8; ++i) {
        v = (v << 8) | static_cast<std::uint64_t>(p[i]);
    }
    return v;
}

std::uint32_t read_uint32_be(const unsigned char* p) {
    return ((static_cast<std::uint32_t>(p[0]) << 24) & UINT32_C(0xFF000000)) |
           ((static_cast<std::uint32_t>(p[1]) << 16) & UINT32_C(0x00FF0000)) |
           ((static_cast<std::uint32_t>(p[2]) <<  8) & UINT32_C(0x0000FF00)) |
           ( static_cast<std::uint32_t>(p[3])        & UINT32_C(0x000000FF));
}

/// Direct chunk read (durin check_for_chunk_read / get_frame_from_chunk):
/// one frame per chunk, native-order unsigned pixels, and either no filter
/// or bitshuffle with LZ4. Anything else returns false for H5Dread.
bool read_frame_chunk(
    hid_t data,
    hsize_t frame_in_block,
    int width,
    int height,
    RawFrame& raw
) {
    const hid_t dcpl = H5Dget_create_plist(data);
    if (dcpl < 0) {
        return false;
    }
    bool ok = false;
    hsize_t cdims[3] = {0, 0, 0};
    if (H5Pget_layout(dcpl) == H5D_CHUNKED &&
        H5Pget_chunk(dcpl, 3, cdims) == 3 && cdims[0] == 1 &&
        cdims[1] == static_cast<hsize_t>(height) &&
        cdims[2] == static_cast<hsize_t>(width)) {
        const int nfilters = H5Pget_nfilters(dcpl);
        bool filter_ok = nfilters == 0;
        unsigned int cd[8] = {0};
        if (nfilters == 1) {
            unsigned int flags = 0;
            std::size_t ncd = 8;
            unsigned int config = 0;
            char fname[32];
            const H5Z_filter_t f = H5Pget_filter2(
                dcpl, 0, &flags, &ncd, cd, sizeof(fname), fname, &config
            );
            // Plain bitshuffle (no LZ4) chunks carry no size header; leave
            // them to the HDF5 filter pipeline.
            filter_ok = f == BSHUF_H5FILTER && ncd > 4 &&
                        cd[4] == BSHUF_H5_COMPRESS_LZ4;
            raw.lz4 = filter_ok;
        }

        const hid_t ftype = H5Dget_type(data);
        const hid_t native = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
        const std::size_t esize = H5Tget_size(ftype);
        const bool type_ok =
            H5Tget_class(ftype) == H5T_INTEGER &&
            H5Tget_sign(ftype) == H5T_SGN_NONE &&
            H5Tget_order(ftype) == H5Tget_order(native) &&
            (esize == 1 || esize == 2 || esize == 4) &&
            (!raw.lz4 || cd[2] == esize);
        H5Tclose(native);
        H5Tclose(ftype);

        const hsize_t offset[3] = {frame_in_block, 0, 0};
        hsize_t nbytes = 0;
        if (filter_ok && type_ok &&
            H5Dget_chunk_storage_size(data, offset, &nbytes) >= 0 &&
            nbytes > 0) {
            raw.chunk.resize(static_cast<std::size_t>(nbytes));
            std::uint32_t skipped = 0;
            if (H5Dread_chunk(data, H5P_DEFAULT, offset, &skipped,
                              raw.chunk.data()) >= 0) {
                // Filter-mask bit 0 set: this chunk was stored unfiltered.
                if ((skipped & UINT32_C(0x00000001)) != 0) {
                    raw.lz4 = false;
                }
                raw.elem_size = esize;
                ok = true;
            }
        }
    }
    H5Pclose(dcpl);
    if (!ok) {
        raw.chunk.clear();
        raw.lz4 = false;
    }
    return ok;
}

/// Locate ``frame_1based`` in the layout and read it. Call with
/// ``hdf5_mutex`` held.
RawFrame read_raw_frame(
    const DataLayout& layout,
    int frame_1based,
    int width,
    int height
) {
    if (frame_1based < 1) {
        throw std::runtime_error("Hdf5ImageReader: frame must be >= 1");
    }
    if (frame_1based > layout.total_frames) {
        throw std::runtime_error(
            "Hdf5ImageReader: frame " + std::to_string(frame_1based) +
            " out of range (nimages=" + std::to_string(layout.total_frames) +
            ")"
        );
    }
    // Block sizes may differ (durin get_dectris_eiger_frame).
    int frame_in_block = frame_1based - 1;
    const DataBlock* block = nullptr;
    for (const auto& b : layout.blocks) {
        if (frame_in_block < b.frames) {
            block = &b;
            break;
        }
        frame_in_block -= b.frames;
    }
    if (block == nullptr) {
        throw std::runtime_error("Hdf5ImageReader: frame not found in data blocks");
    }

    hid_t data = H5Dopen2(layout.group, block->name.c_str(), H5P_DEFAULT);
    if (data < 0) {
        throw std::runtime_error(
            "Hdf5ImageReader: failed to open " + block->name
        );
    }
    RawFrame raw;
    raw.npix =
        static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    if (read_frame_chunk(
            data, static_cast<hsize_t>(frame_in_block), width, height, raw
        )) {
        H5Dclose(data);
        return raw;
    }

    const hid_t ftype = H5Dget_type(data);
    raw.elem_size = H5Tget_size(ftype);
    H5Tclose(ftype);

    hid_t dataspace = H5Dget_space(data);
    const hsize_t offset_in[3] = {
        static_cast<hsize_t>(frame_in_block), 0, 0
    };
    const hsize_t count[3] = {
        1, static_cast<hsize_t>(height), static_cast<hsize_t>(width)
    };
    hid_t memspace = H5Screate_simple(3, count, nullptr);
    if (H5Sselect_hyperslab(
            dataspace, H5S_SELECT_SET, offset_in, nullptr, count, nullptr
        ) < 0) {
        H5Sclose(memspace);
        H5Sclose(dataspace);
        H5Dclose(data);
        throw std::runtime_error("Hdf5ImageReader: hyperslab select failed");
    }
    raw.values.resize(raw.npix);
    if (H5Dread(
            data, H5T_NATIVE_UINT, memspace, dataspace, H5P_DEFAULT,
            raw.values.data()
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
    return raw;
}

/// Unsigned pixel values from a raw chunk (durin bslz4_decompress).
void decode_chunk(RawFrame& raw) {
    const std::size_t nbytes = raw.npix * raw.elem_size;
    std::vector<unsigned char> decoded;
    const unsigned char* src = raw.chunk.data();
    if (raw.lz4) {
        // Header: uncompressed size (uint64 BE), block size in bytes (uint32 BE).
        if (raw.chunk.size() < 12 || read_uint64_be(src) != nbytes) {
            throw std::runtime_error(
                "Hdf5ImageReader: bitshuffle-LZ4 chunk size mismatch"
            );
        }
        const std::size_t block_size = read_uint32_be(src + 8) / raw.elem_size;
        if (block_size == 0) {
            throw std::runtime_error(
                "Hdf5ImageReader: bitshuffle-LZ4 block size is 0"
            );
        }
        decoded.resize(nbytes);
        if (bshuf_decompress_lz4(
                raw.chunk.data() + 12, decoded.data(), raw.npix, raw.elem_size,
                block_size
            ) < 0) {
            throw std::runtime_error(
                "Hdf5ImageReader: bitshuffle-LZ4 decompression failed"
            );
        }
        src = decoded.data();
    } else if (raw.chunk.size() != nbytes) {
        throw std::runtime_error("Hdf5ImageReader: unfiltered chunk size mismatch");
    }

    raw.values.resize(raw.npix);
    if (raw.elem_size == 1) {
        for (std::size_t i = 0; i < raw.npix; ++i) {
            raw.values[i] = src[i];
        }
    } else if (raw.elem_size == 2) {
        for (std::size_t i = 0; i < raw.npix; ++i) {
            std::uint16_t v;
            std::memcpy(&v, src + 2 * i, 2);
            raw.values[i] = v;
        }
    } else {
        std::memcpy(raw.values.data(), src, nbytes);
    }
    raw.chunk.clear();
    raw.chunk.shrink_to_fit();
}

/// Decode ``raw`` into MOSFLM pixels. Runs without ``hdf5_mutex``.
void finish_pixels(
    RawFrame& raw,
    unsigned int error_val,
    const PixelMask* mask,
    std::vector<std::int32_t>& pixels
) {
    if (!raw.chunk.empty()) {
        decode_chunk(raw);
    }
    // All bits set at the stored width is DECTRIS's "no valid value"; it
    // means masked only when it lies above the saturation level (otherwise
    // a gap and a saturated pixel are indistinguishable, and DEFPIX masks
    // pixels that are overloaded throughout).
    const unsigned int all_ones =
        raw.elem_size == 1 ? UINT32_C(0x000000FF)
        : raw.elem_size == 2 ? UINT32_C(0x0000FFFF)
                             : UINT32_C(0xFFFFFFFF);
    const bool all_ones_masked =
        all_ones == UINT32_C(0xFFFFFFFF) || error_val < all_ones;

    pixels.resize(raw.npix);
    // Sentinels: −1 = masked (pixel_mask, module gap / bad pixel);
    // −2 = overloaded (≥ saturation value). Integrators skip both. durin
    // gives some mask bits −2 (as Neggia); here every mask bit means −1 so
    // a bad pixel is never reported as an overload.
    for (std::size_t i = 0; i < raw.npix; ++i) {
        const unsigned int v = raw.values[i];
        if ((mask != nullptr && (*mask)[i] != 0) ||
            (all_ones_masked && v == all_ones)) {
            pixels[i] = -1;
        } else if (v >= error_val) {
            pixels[i] = -2;
        } else {
            pixels[i] = static_cast<std::int32_t>(v);
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
    DiffractionImage out;
    RawFrame raw;
    std::shared_ptr<const PixelMask> mask;
    unsigned int error_val = 0;
    {
        const std::lock_guard<std::mutex> lock(hdf5_mutex());
        const hid_t hdf = open_file(path);
        DataLayout layout;
        try {
            layout = open_data_layout(hdf);
            std::string det = find_nx_group(
                hdf, "/entry/instrument", "/entry/instrument/detector",
                "NXdetector"
            );
            if (det.empty()) {
                det = "/entry/instrument/detector";
            }
            MasterMeta meta = read_metadata(hdf, det, frame, layout);
            if (meta.nimages <= 0) {
                meta.nimages = layout.total_frames;
            }

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
            out.metadata.saturation_counts = meta.error_val;
            error_val = meta.error_val;

            mask = read_pixel_mask(hdf, path, det, meta.width, meta.height);
            raw = read_raw_frame(layout, frame, meta.width, meta.height);

            H5Gclose(layout.group);
            H5Fclose(hdf);
        } catch (...) {
            if (layout.group >= 0) {
                H5Gclose(layout.group);
            }
            H5Fclose(hdf);
            throw;
        }
    }
    finish_pixels(raw, error_val, mask.get(), out.pixels);
    return out;
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
