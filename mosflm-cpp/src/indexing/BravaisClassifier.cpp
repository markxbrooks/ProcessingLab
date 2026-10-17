#include "mosflm/indexing/BravaisClassifier.hpp"

#include <cmath>
#include <cstring>
#include <stdexcept>

// DPS / Cambridge headers. util.h must precede bravais.h so matrix_3x3 exists
// before bravais_goof.h. Linked from lib/dps_index.a + lib/libdps.a.
extern "C" {
#include "util.h"
#include "bravais.h"
}

namespace mosflm::indexing {
namespace {

struct matrix_3x3 eigen_to_dps(const Eigen::Matrix3d& M) {
    struct matrix_3x3 om{};
    for (int i = 0; i < 3; ++i) {
        om.r[i].x = static_cast<float>(M(i, 0));
        om.r[i].y = static_cast<float>(M(i, 1));
        om.r[i].z = static_cast<float>(M(i, 2));
    }
    return om;
}

Eigen::Matrix3d dps_to_eigen(const struct matrix_3x3& om) {
    Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
    for (int i = 0; i < 3; ++i) {
        M(i, 0) = static_cast<double>(om.r[i].x);
        M(i, 1) = static_cast<double>(om.r[i].y);
        M(i, 2) = static_cast<double>(om.r[i].z);
    }
    return M;
}

}  // namespace

std::vector<BravaisCandidate> BravaisClassifier::classify(
    const Eigen::Matrix3d& reciprocal_om,
    int output_type
) const {
    struct matrix_3x3 om = eigen_to_dps(reciprocal_om);
    struct matrix_3x3 cl_om[45];
    float crit[45];
    char sys_id[45][3];
    std::memset(cl_om, 0, sizeof(cl_om));
    std::memset(crit, 0, sizeof(crit));
    std::memset(sys_id, 0, sizeof(sys_id));

    ::bravais(&om, cl_om, crit, sys_id, output_type);

    std::vector<BravaisCandidate> out;
    out.reserve(44);
    for (int i = 1; i <= 44; ++i) {
        BravaisCandidate c;
        c.character = i;
        c.system_id.assign(sys_id[i], sys_id[i] + 2);
        c.criterion = crit[i];
        c.orientation = dps_to_eigen(cl_om[i]);
        out.push_back(std::move(c));
    }
    return out;
}

std::vector<BravaisCandidate> BravaisClassifier::classify_orthogonal_cell(
    double a_A,
    double b_A,
    double c_A
) const {
    if (!(a_A > 0.0 && b_A > 0.0 && c_A > 0.0)) {
        throw std::invalid_argument(
            "BravaisClassifier::classify_orthogonal_cell: non-positive edge"
        );
    }
    Eigen::Matrix3d om = Eigen::Matrix3d::Zero();
    om(0, 0) = 1.0 / a_A;
    om(1, 1) = 1.0 / b_A;
    om(2, 2) = 1.0 / c_A;
    return classify(om, 1);
}

}  // namespace mosflm::indexing
