#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/crystal/SpaceGroupInfo.hpp"
#include "test_util.hpp"

using mosflm::UnitCell;
using mosflm::crystal::CrystalMetadata;
using mosflm::crystal::SpaceGroupInfo;

int main() {
    const CrystalMetadata meta =
        CrystalMetadata::from_parameters(58.0, 58.0, 150.0, 90.0, 90.0, 90.0,
                                         "P 4");
    CHECK(meta.spacegroup != nullptr);
    CHECK_EQ(SpaceGroupInfo::number(meta.spacegroup), 75);
    CHECK(!SpaceGroupInfo::hall(meta.spacegroup).empty());

    // Roundtrip to the MOSFLM cell type.
    const UnitCell cell = meta.to_mosflm_cell();
    CHECK_CLOSE(cell.a, 58.0, 1e-9);
    CHECK_CLOSE(cell.c, 150.0, 1e-9);
    CHECK_CLOSE(cell.gamma_deg, 90.0, 1e-9);

    // Convenience constructor from the MOSFLM cell.
    const CrystalMetadata meta2 = CrystalMetadata::from_mosflm_cell(cell, "P 1");
    CHECK(meta2.spacegroup != nullptr);
    CHECK_EQ(SpaceGroupInfo::number(meta2.spacegroup), 1);

    // Reciprocal cell from Gemmi.
    const auto recip = meta.reciprocal();
    CHECK_CLOSE(recip.a, 1.0 / 58.0, 1e-9);

    // describe() is non-empty and reports the space group.
    const std::string text = meta.describe();
    CHECK(!text.empty());
    CHECK(text.find("a=58") != std::string::npos);

    // Space-group lookup helpers.
    const auto sg = SpaceGroupInfo::find_by_name("P 21 21 21");
    CHECK(sg != nullptr);
    CHECK_EQ(SpaceGroupInfo::number(sg), 19);

    CHECK(SpaceGroupInfo::find_by_name("definitely not a space group") ==
          nullptr);
    CHECK_EQ(SpaceGroupInfo::number(nullptr), 0);
    CHECK(SpaceGroupInfo::hall(nullptr).empty());

    // Default space group is P 1.
    const CrystalMetadata p1 =
        CrystalMetadata::from_parameters(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    CHECK_EQ(SpaceGroupInfo::number(p1.spacegroup), 1);

    return mosflm_test::finish("gemmi_crystal");
}
