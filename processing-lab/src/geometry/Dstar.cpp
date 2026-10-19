#include "mosflm/geometry/Dstar.hpp"

#include <cmath>
#include <stdexcept>

namespace mosflm {

Eigen::Vector3d Dstar::q_vector(
    const OrientationMatrix& ub,
    int h,
    int k,
    int l
) {
    return ub.miller_to_lab(h, k, l);
}

double Dstar::d_spacing(const Eigen::Vector3d& q) {
    const double len = q.norm();
    if (len <= 0.0) {
        return 0.0;
    }
    return 1.0 / len;
}

double Dstar::d_spacing(const UnitCell& cell, int h, int k, int l) {
    const OrientationMatrix ub = OrientationMatrix::from_cell_and_rotation(cell);
    return d_spacing(q_vector(ub, h, k, l));
}

double Dstar::excitation_error(
    const Eigen::Vector3d& q,
    double wavelength_A
) {
    if (!(wavelength_A > 0.0)) {
        throw std::invalid_argument("Dstar::excitation_error: invalid wavelength");
    }
    const double inv_lambda = 1.0 / wavelength_A;
    const Eigen::Vector3d s0(0.0, 0.0, inv_lambda);
    return (s0 + q).norm() - inv_lambda;
}

}  // namespace mosflm
