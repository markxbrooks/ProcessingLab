#include "mosflm/geometry/Rotmat.hpp"

#include <cmath>

namespace mosflm {
namespace {

// MOSFLM trigparam.f DTOR / RTOD (Clemens Vonrhein constants).
constexpr double kDtor = 0.0174532925;
constexpr double kRtod = 57.2957795131;

}  // namespace

Eigen::Matrix3d rotmat_from_missetting(
    double phix,
    double phiy,
    double phiz,
    bool degrees
) {
    double ax = phix;
    double ay = phiy;
    double az = phiz;
    if (degrees) {
        ax *= kDtor;
        ay *= kDtor;
        az *= kDtor;
    }
    const double cx = std::cos(ax);
    const double sx = std::sin(ax);
    const double cy = std::cos(ay);
    const double sy = std::sin(ay);
    const double cz = std::cos(az);
    const double sz = std::sin(az);

    // [R] = [phiz] . [phiy] . [phix]  (column-vector convention, R * v)
    Eigen::Matrix3d Rx = Eigen::Matrix3d::Identity();
    Rx(1, 1) = cx;
    Rx(1, 2) = -sx;
    Rx(2, 1) = sx;
    Rx(2, 2) = cx;

    Eigen::Matrix3d Ry = Eigen::Matrix3d::Identity();
    Ry(0, 0) = cy;
    Ry(0, 2) = sy;
    Ry(2, 0) = -sy;
    Ry(2, 2) = cy;

    Eigen::Matrix3d Rz = Eigen::Matrix3d::Identity();
    Rz(0, 0) = cz;
    Rz(0, 1) = -sz;
    Rz(1, 0) = sz;
    Rz(1, 1) = cz;

    (void)kRtod;
    return Rz * Ry * Rx;
}

}  // namespace mosflm
