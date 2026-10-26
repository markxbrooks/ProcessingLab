// Hdf5ImageReader on small synthetic files: DECTRIS data_NNNNNN blocks of
// unequal size, the NXmx single ``data`` dataset (plain and virtual), an
// NXdetector found by NX_class, ``units`` scaling, ``pixel_mask``, and the
// direct bitshuffle-LZ4 chunk read beside the HDF5 filter pipeline.
#include "mosflm/image/DiffractionImage.hpp"
#include "test_util.hpp"

#include <hdf5.h>
#include <hdf5_hl.h>

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

extern "C" {
#include "bshuf_h5filter.h"
}

namespace {

namespace fs = std::filesystem;
using mosflm::image::Hdf5ImageReader;

constexpr int kW = 6;
constexpr int kH = 4;

enum class Storage { BitshuffleLz4, Bitshuffle, Chunked, Contiguous };

/// Pixel value written for (frame, index); frame is 1-based over the sweep.
std::uint32_t value_at(int frame, int i) {
    return static_cast<std::uint32_t>(frame * 100 + i);
}

hid_t make_groups(hid_t file, const std::string& path, const char* nx_class) {
    const hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(lcpl, 1);
    const hid_t g = H5Gcreate2(file, path.c_str(), lcpl, H5P_DEFAULT, H5P_DEFAULT);
    H5Pclose(lcpl);
    H5LTset_attribute_string(file, path.c_str(), "NX_class", nx_class);
    return g;
}

void write_scalar(hid_t loc, const char* name, double v, const char* units,
                  bool vlen_units = false) {
    const hsize_t one = 1;
    H5LTmake_dataset_double(loc, name, 0, &one, &v);
    if (units == nullptr) {
        return;
    }
    if (!vlen_units) {
        H5LTset_attribute_string(loc, name, "units", units);
        return;
    }
    const hid_t ds = H5Dopen2(loc, name, H5P_DEFAULT);
    const hid_t type = H5Tcopy(H5T_C_S1);
    H5Tset_size(type, H5T_VARIABLE);
    const hid_t space = H5Screate(H5S_SCALAR);
    const hid_t attr = H5Acreate2(ds, "units", type, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, type, &units);
    H5Aclose(attr);
    H5Sclose(space);
    H5Tclose(type);
    H5Dclose(ds);
}

void write_int(hid_t loc, const char* name, int v) {
    H5LTmake_dataset_int(loc, name, 0, nullptr, &v);
}

void write_mask(hid_t loc, const char* name, const std::vector<std::uint32_t>& m) {
    const hsize_t dims[2] = {kH, kW};
    const hid_t space = H5Screate_simple(2, dims, nullptr);
    const hid_t ds = H5Dcreate2(loc, name, H5T_STD_U32LE, space, H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(ds, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, m.data());
    H5Dclose(ds);
    H5Sclose(space);
}

/// (n, kH, kW) image dataset holding sweep frames first .. first+n−1, with
/// ``patch`` overriding pixel values (index → value) in every frame.
void write_frames(hid_t loc, const char* name, hid_t file_type, int first, int n,
                  Storage storage,
                  const std::vector<std::pair<int, std::uint32_t>>& patch = {}) {
    const hsize_t dims[3] = {static_cast<hsize_t>(n), kH, kW};
    const hid_t space = H5Screate_simple(3, dims, nullptr);
    const hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    if (storage != Storage::Contiguous) {
        const hsize_t chunk[3] = {1, kH, kW};
        H5Pset_chunk(dcpl, 3, chunk);
    }
    if (storage == Storage::BitshuffleLz4 || storage == Storage::Bitshuffle) {
        // {block size (0 = default), compression}; set_local prepends the
        // version and element size.
        const unsigned int cd[2] = {
            0, storage == Storage::BitshuffleLz4 ? BSHUF_H5_COMPRESS_LZ4 : 0u};
        H5Pset_filter(dcpl, BSHUF_H5FILTER, H5Z_FLAG_MANDATORY, 2, cd);
    }
    const hid_t ds = H5Dcreate2(loc, name, file_type, space, H5P_DEFAULT, dcpl,
                                H5P_DEFAULT);
    std::vector<std::uint32_t> data(static_cast<std::size_t>(n) * kW * kH);
    for (int f = 0; f < n; ++f) {
        for (int i = 0; i < kW * kH; ++i) {
            data[static_cast<std::size_t>(f) * kW * kH + i] = value_at(first + f, i);
        }
        for (const auto& [i, v] : patch) {
            data[static_cast<std::size_t>(f) * kW * kH + i] = v;
        }
    }
    H5Dwrite(ds, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    H5Dclose(ds);
    H5Pclose(dcpl);
    H5Sclose(space);
}

/// DECTRIS layout: data_000001 (3 frames) + data_000002 (2 frames), uint16,
/// bitshuffle-LZ4 (or another storage), units in µm / cm, mask under
/// detectorSpecific, saturation 60000.
fs::path write_dectris(const fs::path& dir, Storage storage) {
    const fs::path path = dir / "dectris_master.h5";
    const hid_t f = H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                              H5P_DEFAULT);
    H5Gclose(make_groups(f, "/entry", "NXentry"));
    const hid_t det = make_groups(f, "/entry/instrument/detector", "NXdetector");
    write_scalar(det, "x_pixel_size", 75.0, "um");
    write_scalar(det, "y_pixel_size", 75.0, "um");
    write_scalar(det, "detector_distance", 25.0, "cm");
    write_int(det, "saturation_value", 60000);
    const hid_t spec = make_groups(f, "/entry/instrument/detector/detectorSpecific",
                                   "NXcollection");
    std::vector<std::uint32_t> mask(kW * kH, 0);
    mask[0] = UINT32_C(0x00000001);   // gap
    mask[8] = UINT32_C(0x00000008);   // hot
    mask[14] = UINT32_C(0x00000100);  // above the low byte: kept
    write_mask(spec, "pixel_mask", mask);
    H5Gclose(spec);
    H5Gclose(det);
    const hid_t data = make_groups(f, "/entry/data", "NXdata");
    // Pixel 20 saturates; 21 is 0xFFFF (masked: above saturation).
    const std::vector<std::pair<int, std::uint32_t>> patch = {
        {20, 60000}, {21, 65535}};
    write_frames(data, "data_000001", H5T_STD_U16LE, 1, 3, storage, patch);
    write_frames(data, "data_000002", H5T_STD_U16LE, 4, 2, storage, patch);
    H5Gclose(data);
    H5Fclose(f);
    return path;
}

/// NXmx layout: detector group "eiger" found by NX_class, one ``data``
/// dataset (or a virtual one over two source datasets in another file),
/// uint32, metres given as a variable-length units string, mask in the
/// detector group.
fs::path write_nxmx(const fs::path& dir, bool virtual_data) {
    const fs::path path = dir / (virtual_data ? "vds.nxs" : "nxmx.nxs");
    if (virtual_data) {
        const fs::path src = dir / "vds_source.h5";
        const hid_t s = H5Fcreate(src.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                                  H5P_DEFAULT);
        write_frames(s, "a", H5T_STD_U32LE, 1, 2, Storage::BitshuffleLz4);
        write_frames(s, "b", H5T_STD_U32LE, 3, 3, Storage::BitshuffleLz4);
        H5Fclose(s);
    }
    const hid_t f = H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                              H5P_DEFAULT);
    H5Gclose(make_groups(f, "/entry", "NXentry"));
    const hid_t det = make_groups(f, "/entry/instrument/eiger", "NXdetector");
    write_scalar(det, "x_pixel_size", 0.000172, "m", true);
    write_scalar(det, "y_pixel_size", 0.000172, "m", true);
    write_scalar(det, "distance", 0.3, "m", true);
    std::vector<std::uint32_t> mask(kW * kH, 0);
    mask[5] = UINT32_C(0x00000002);  // dead
    write_mask(det, "pixel_mask", mask);
    H5Gclose(det);
    const hid_t data = make_groups(f, "/entry/data", "NXdata");
    if (!virtual_data) {
        write_frames(data, "data", H5T_STD_U32LE, 1, 5, Storage::Chunked);
    } else {
        const hsize_t dims[3] = {5, kH, kW};
        const hid_t vspace = H5Screate_simple(3, dims, nullptr);
        const hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        const char* names[2] = {"a", "b"};
        const hsize_t counts[2] = {2, 3};
        hsize_t start = 0;
        for (int k = 0; k < 2; ++k) {
            const hsize_t sdims[3] = {counts[k], kH, kW};
            const hid_t sspace = H5Screate_simple(3, sdims, nullptr);
            const hsize_t off[3] = {start, 0, 0};
            H5Sselect_hyperslab(vspace, H5S_SELECT_SET, off, nullptr, sdims, nullptr);
            H5Pset_virtual(dcpl, vspace, "vds_source.h5", names[k], sspace);
            H5Sclose(sspace);
            start += counts[k];
        }
        H5Sselect_all(vspace);
        const hid_t ds = H5Dcreate2(data, "data", H5T_STD_U32LE, vspace,
                                    H5P_DEFAULT, dcpl, H5P_DEFAULT);
        H5Dclose(ds);
        H5Pclose(dcpl);
        H5Sclose(vspace);
    }
    H5Gclose(data);
    H5Fclose(f);
    return path;
}

void check_frame(const fs::path& path, int frame, const std::string& label,
                 const std::vector<int>& masked, const std::vector<int>& overloaded) {
    const auto img = Hdf5ImageReader{}.read(path, frame);
    CHECK_MSG(img.metadata.width == kW && img.metadata.height == kH, label + " size");
    CHECK_MSG(img.pixels.size() == static_cast<std::size_t>(kW * kH), label + " pixels");
    if (img.pixels.size() != static_cast<std::size_t>(kW * kH)) {
        return;
    }
    for (int i = 0; i < kW * kH; ++i) {
        std::int32_t want = static_cast<std::int32_t>(value_at(frame, i));
        for (int m : masked) {
            if (m == i) want = -1;
        }
        for (int o : overloaded) {
            if (o == i) want = -2;
        }
        CHECK_MSG(img.pixels[static_cast<std::size_t>(i)] == want,
                  label + " frame " + std::to_string(frame) + " pixel " +
                      std::to_string(i) + ": got " +
                      std::to_string(img.pixels[static_cast<std::size_t>(i)]) +
                      " want " + std::to_string(want));
    }
}

}  // namespace

int main() {
    bshuf_register_h5filter();
    const fs::path dir = mosflm_test::output_dir("hdf5_reader");

    // DECTRIS blocks: every storage gives the same pixels. Mask bits 0 and 3
    // mask; bit 8 does not. 0xFFFF lies above saturation 60000, so it is a
    // masked sentinel; 60000 itself is an overload.
    const std::pair<Storage, const char*> storages[] = {
        {Storage::BitshuffleLz4, "bslz4"},
        {Storage::Bitshuffle, "bitshuffle"},
        {Storage::Chunked, "chunked"},
        {Storage::Contiguous, "contiguous"},
    };
    for (const auto& [storage, label] : storages) {
        const fs::path sub = dir / label;
        fs::create_directories(sub);
        const fs::path p = write_dectris(sub, storage);
        CHECK_MSG(Hdf5ImageReader{}.num_images(p) == 5,
                  std::string(label) + " num_images");
        for (int frame : {1, 3, 4, 5}) {
            check_frame(p, frame, label, {0, 8, 21}, {20});
        }
        const auto img = Hdf5ImageReader{}.read(p, 2);
        CHECK_CLOSE(img.metadata.pixel_size_x_mm, 0.075, 1e-12);
        CHECK_CLOSE(img.metadata.detector_distance_mm, 250.0, 1e-9);
        CHECK(img.metadata.saturation_counts == 60000u);
        bool threw = false;
        try {
            (void)Hdf5ImageReader{}.read(p, 6);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        CHECK_MSG(threw, std::string(label) + " frame past the end throws");
    }

    // NXmx single dataset and its virtual-dataset form.
    for (const bool vds : {false, true}) {
        const fs::path sub = dir / (vds ? "vds" : "nxmx");
        fs::create_directories(sub);
        const fs::path p = write_nxmx(sub, vds);
        const std::string label = vds ? "vds" : "nxmx";
        CHECK_MSG(Hdf5ImageReader{}.num_images(p) == 5, label + " num_images");
        for (int frame : {1, 2, 3, 5}) {
            check_frame(p, frame, label, {5}, {});
        }
        const auto img = Hdf5ImageReader{}.read(p, 1);
        CHECK_CLOSE(img.metadata.pixel_size_x_mm, 0.172, 1e-12);
        CHECK_CLOSE(img.metadata.detector_distance_mm, 300.0, 1e-9);
    }

    return mosflm_test::finish("hdf5_reader");
}
