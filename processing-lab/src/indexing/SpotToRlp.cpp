#include "mosflm/indexing/SpotToRlp.hpp"

#include <cmath>

namespace mosflm::indexing {

std::size_t spots_to_rlp_xyz(
    const SpotList& spots,
    const std::vector<std::size_t>& keep,
    const SpotToRlpGeometry& geometry,
    std::vector<float>& x,
    std::vector<float>& y,
    std::vector<float>& z
) {
    x.clear();
    y.clear();
    z.clear();
    if (!(geometry.distance_mm > 0.0) || !(geometry.wavelength_A > 0.0)) {
        return 0;
    }

    std::vector<std::size_t> idx = keep;
    if (idx.empty()) {
        idx.resize(spots.spots.size());
        for (std::size_t i = 0; i < idx.size(); ++i) {
            idx[i] = i;
        }
    }

    constexpr double kDeg2Rad = 3.14159265358979323846 / 180.0;
    const double dist = geometry.distance_mm;
    const double inv_l = 1.0 / geometry.wavelength_A;
    const double omega = spots.header.omega_deg * kDeg2Rad;
    const double finvertx = (spots.header.invflag != 0) ? -1.0 : 1.0;
    const double fr11 = std::cos(omega) * finvertx;
    const double fr12 = std::sin(omega);
    const double fr21 = -std::sin(omega) * finvertx;
    const double fr22 = std::cos(omega);
    const double beam_x = geometry.beam_x_mm;
    const double beam_y = geometry.beam_y_mm;

    x.reserve(idx.size());
    y.reserve(idx.size());
    z.reserve(idx.size());
    for (std::size_t ki : idx) {
        if (ki >= spots.spots.size()) {
            continue;
        }
        const auto& s = spots.spots[ki];
        const double xt = s.x_mm - beam_x;
        const double yt = s.y_mm - beam_y;
        const double xs = xt * fr11 + yt * fr12;
        const double ys = xt * fr21 + yt * fr22;
        const double rlen = std::sqrt(xs * xs + ys * ys + dist * dist);
        if (!(rlen > 1e-12)) {
            continue;
        }
        const double f = inv_l / rlen;
        const double qx = f * xs;
        const double qy = f * ys;
        const double qz = f * dist - inv_l;
        const double phi_rad = s.phi_end_deg * kDeg2Rad;
        const double c = std::cos(phi_rad);
        const double sn = std::sin(phi_rad);
        // ROTMATY (TRANSP∘RTMATS about Y): x'=c x − s z, z'=s x + c z
        x.push_back(static_cast<float>(c * qx - sn * qz));
        y.push_back(static_cast<float>(qy));
        z.push_back(static_cast<float>(sn * qx + c * qz));
    }
    return x.size();
}

}  // namespace mosflm::indexing
