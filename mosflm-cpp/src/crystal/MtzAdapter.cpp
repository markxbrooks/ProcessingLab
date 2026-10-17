#include "mosflm/crystal/MtzAdapter.hpp"

#include <gemmi/mtz.hpp>

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

}  // namespace mosflm::crystal
