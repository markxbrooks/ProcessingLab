#include "mosflm/geometry/LabFrame.hpp"
#include "test_util.hpp"

int main() {
    const Eigen::Matrix3d R = mosflm::mosflm_lab_to_cpp();
    CHECK_CLOSE(R.determinant(), 1.0, 1e-12);

    // Photon travel: MOSFLM +X → C++ +Z
    const Eigen::Vector3d beam_dir_m(1.0, 0.0, 0.0);
    const Eigen::Vector3d beam_dir_c = R * beam_dir_m;
    CHECK_CLOSE(beam_dir_c.x(), 0.0, 1e-12);
    CHECK_CLOSE(beam_dir_c.y(), 0.0, 1e-12);
    CHECK_CLOSE(beam_dir_c.z(), 1.0, 1e-12);

    // MOSFLM S0 (antiparallel to beam) → C++ −Z
    const Eigen::Vector3d s0_m(-1.0, 0.0, 0.0);
    const Eigen::Vector3d s0_c = R * s0_m;
    CHECK_CLOSE(s0_c.x(), 0.0, 1e-12);
    CHECK_CLOSE(s0_c.y(), 0.0, 1e-12);
    CHECK_CLOSE(s0_c.z(), -1.0, 1e-12);

    // Spindle: MOSFLM +Z → C++ +X
    const Eigen::Vector3d spin_m(0.0, 0.0, 1.0);
    const Eigen::Vector3d spin_c = R * spin_m;
    CHECK_CLOSE(spin_c.x(), 1.0, 1e-12);
    CHECK_CLOSE(spin_c.y(), 0.0, 1e-12);
    CHECK_CLOSE(spin_c.z(), 0.0, 1e-12);

    // Identity A → same R
    const Eigen::Matrix3d Ub = mosflm::mosflm_amat_to_cpp_ub(
        Eigen::Matrix3d::Identity()
    );
    CHECK_CLOSE((Ub - R).norm(), 0.0, 1e-12);

    // Detector: XDET↔y_mm, YDET↔x_mm relative to BEAM
    double xdet = 0.0;
    double ydet = 0.0;
    mosflm::cpp_mm_to_mosflm_xdet(
        10.0, -5.0, 94.0, 94.7, 0.0816, 0.0816, xdet, ydet
    );
    CHECK_CLOSE(xdet, 94.0 / 0.0816 - 5.0 / 0.0816, 1e-9);
    CHECK_CLOSE(ydet, 94.7 / 0.0816 + 10.0 / 0.0816, 1e-9);

    return mosflm_test::finish("lab_frame");
}
