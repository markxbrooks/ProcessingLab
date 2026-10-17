#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/geometry/UnitCell.hpp"
#include "mosflm/legacy/LegacyMosflmAdapter.hpp"

#include <iostream>

int main() {
    std::cout << "mosflm-cpp " << "0.1.0" << '\n';
    std::cout << mosflm::legacy::LegacyMosflmAdapter::oracle_hint() << '\n';

    mosflm::UnitCell cell{58.0, 58.0, 150.0, 90.0, 90.0, 90.0};
    std::cout << "demo cell: " << cell.describe() << '\n';

    const auto crystal =
        mosflm::crystal::CrystalMetadata::from_mosflm_cell(cell, "P 41 21 2");
    std::cout << "gemmi crystal: " << crystal.describe() << '\n';

    const auto ub = mosflm::OrientationMatrix::from_cell_and_rotation(cell);
    const Eigen::Vector3d r = ub.miller_to_lab(1, 0, 0);
    std::cout << "UB*(1,0,0) = [" << r.transpose() << "] Å^-1\n";

    mosflm::DetectorGeometry det;
    det.distance_mm = 250.0;
    std::cout << "detector distance " << det.distance_mm << " mm\n";
    return 0;
}
