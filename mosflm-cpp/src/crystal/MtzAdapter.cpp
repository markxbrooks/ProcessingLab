#include "mosflm/crystal/MtzAdapter.hpp"

#include <gemmi/mtz.hpp>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace mosflm::crystal {

void MtzAdapter::write_simple(
    const std::filesystem::path& path,
    const CrystalMetadata& crystal,
    const std::vector<SimpleReflection>& reflections
) {
    gemmi::Mtz mtz;
    mtz.title = "mosflm-cpp simple";
    mtz.cell = crystal.cell;
    if (crystal.spacegroup != nullptr) {
        mtz.spacegroup = crystal.spacegroup;
        mtz.spacegroup_number = crystal.spacegroup->number;
        mtz.spacegroup_name = crystal.spacegroup->hm;
    }

    mtz.add_base();  // H K L on dataset 0
    const gemmi::Mtz::Dataset& ds = mtz.add_dataset("mosflm");
    mtz.add_column("I", 'J', ds.id, -1, false);
    mtz.add_column("SIGI", 'Q', ds.id, -1, false);

    std::vector<float> data;
    data.reserve(reflections.size() * 5);
    for (const auto& r : reflections) {
        data.push_back(static_cast<float>(r.h));
        data.push_back(static_cast<float>(r.k));
        data.push_back(static_cast<float>(r.l));
        data.push_back(r.intensity);
        data.push_back(r.sigma);
    }
    mtz.set_data(data.data(), data.size());

    try {
        mtz.write_to_file(path.string());
    } catch (const std::exception& ex) {
        throw std::runtime_error(
            std::string("MtzAdapter::write_simple failed: ") + ex.what()
        );
    }
}

std::vector<SimpleReflection> MtzAdapter::read_simple(
    const std::filesystem::path& path
) {
    gemmi::Mtz mtz;
    try {
        mtz.read_file(path.string());
    } catch (const std::exception& ex) {
        throw std::runtime_error(
            std::string("MtzAdapter::read_simple failed: ") + ex.what()
        );
    }
    const gemmi::Mtz::Column* col_h = mtz.column_with_label("H");
    const gemmi::Mtz::Column* col_k = mtz.column_with_label("K");
    const gemmi::Mtz::Column* col_l = mtz.column_with_label("L");
    const gemmi::Mtz::Column* col_i = mtz.column_with_label("I");
    const gemmi::Mtz::Column* col_s = mtz.column_with_label("SIGI");
    if (!col_h || !col_k || !col_l || !col_i || !col_s) {
        throw std::runtime_error(
            "MtzAdapter::read_simple: missing H/K/L/I/SIGI columns"
        );
    }
    std::vector<SimpleReflection> out;
    out.reserve(static_cast<std::size_t>(mtz.nreflections));
    for (int i = 0; i < mtz.nreflections; ++i) {
        SimpleReflection r;
        r.h = static_cast<int>(std::lround((*col_h)[static_cast<std::size_t>(i)]));
        r.k = static_cast<int>(std::lround((*col_k)[static_cast<std::size_t>(i)]));
        r.l = static_cast<int>(std::lround((*col_l)[static_cast<std::size_t>(i)]));
        r.intensity = (*col_i)[static_cast<std::size_t>(i)];
        r.sigma = (*col_s)[static_cast<std::size_t>(i)];
        out.push_back(r);
    }
    return out;
}

}  // namespace mosflm::crystal
