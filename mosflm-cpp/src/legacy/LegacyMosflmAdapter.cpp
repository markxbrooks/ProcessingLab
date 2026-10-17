#include "mosflm/legacy/LegacyMosflmAdapter.hpp"

namespace mosflm::legacy {

void LegacyMosflmAdapter::to_cell6(const UnitCell& cell, float out[6]) {
    out[0] = static_cast<float>(cell.a);
    out[1] = static_cast<float>(cell.b);
    out[2] = static_cast<float>(cell.c);
    out[3] = static_cast<float>(cell.alpha_deg);
    out[4] = static_cast<float>(cell.beta_deg);
    out[5] = static_cast<float>(cell.gamma_deg);
}

UnitCell LegacyMosflmAdapter::from_cell6(const float in[6]) {
    UnitCell cell;
    cell.a = in[0];
    cell.b = in[1];
    cell.c = in[2];
    cell.alpha_deg = in[3];
    cell.beta_deg = in[4];
    cell.gamma_deg = in[5];
    return cell;
}

std::string LegacyMosflmAdapter::oracle_hint() {
    return "Run ../bin/ipmosflm with DYLD_LIBRARY_PATH=$CCP4/lib:"
           "/opt/homebrew/opt/hdf5/lib (see porting.md).";
}

}  // namespace mosflm::legacy
