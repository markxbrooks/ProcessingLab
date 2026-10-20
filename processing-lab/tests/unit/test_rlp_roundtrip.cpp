/**
 * Spot → RLP round-trip: MOSFLM FR_MAT / triangular / ROTMATY recovers A·hkl.
 */

#include "mosflm/geometry/LabFrame.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/geometry/UnitCell.hpp"
#include "mosflm/indexing/SpotList.hpp"
#include "mosflm/indexing/SpotToRlp.hpp"
#include "test_util.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <random>
#include <vector>

using mosflm::OrientationMatrix;
using mosflm::UnitCell;
using mosflm::cpp_ub_to_mosflm_amat;
using mosflm::indexing::Spot;
using mosflm::indexing::SpotList;
using mosflm::indexing::SpotToRlpGeometry;
using mosflm::indexing::spots_to_rlp_xyz;

namespace {

constexpr double kDeg2Rad = 3.14159265358979323846 / 180.0;

/// q = ROTMATY⁻¹(φ)·r0 must lie on the Ewald sphere: |q+(0,0,1/λ)| = 1/λ.
bool on_ewald(const Eigen::Vector3d& q, double inv_l, double tol = 1e-4) {
    const Eigen::Vector3d S(q.x(), q.y(), q.z() + inv_l);
    return std::fabs(S.norm() - inv_l) <= tol * inv_l;
}

bool rlp_to_spot_mm(
    const Eigen::Vector3d& r0,
    double phi_deg,
    double distance_mm,
    double wavelength_A,
    double beam_x_mm,
    double beam_y_mm,
    double& x_mm,
    double& y_mm
) {
    const double inv_l = 1.0 / wavelength_A;
    const double c = std::cos(phi_deg * kDeg2Rad);
    const double sn = std::sin(phi_deg * kDeg2Rad);
    // Inverse ROTMATY so spots_to_rlp recovers r0.
    const double qx = c * r0.x() + sn * r0.z();
    const double qy = r0.y();
    const double qz = -sn * r0.x() + c * r0.z();
    if (!on_ewald(Eigen::Vector3d(qx, qy, qz), inv_l)) {
        return false;
    }
    const double f = (qz + inv_l) / distance_mm;
    if (!(f > 1e-12)) {
        return false;
    }
    x_mm = qx / f + beam_x_mm;
    y_mm = qy / f + beam_y_mm;
    return std::isfinite(x_mm) && std::isfinite(y_mm);
}

}  // namespace

int main() {
    UnitCell cell;
    cell.a = 58.0;
    cell.b = 58.0;
    cell.c = 150.0;
    const OrientationMatrix orient =
        OrientationMatrix::from_cell_and_rotation(cell);
    const double wave = 1.0;
    const Eigen::Matrix3d amat = cpp_ub_to_mosflm_amat(orient.ub(), wave);

    SpotToRlpGeometry geo;
    geo.distance_mm = 200.0;
    geo.wavelength_A = wave;
    geo.beam_x_mm = 160.0;
    geo.beam_y_mm = 160.0;

    SpotList spots;
    spots.header.omega_deg = 0.0;
    spots.header.invflag = 0;
    spots.header.beam_x_mm = geo.beam_x_mm;
    spots.header.beam_y_mm = geo.beam_y_mm;
    spots.header.pixel_mm = 0.075;

    std::vector<Eigen::Vector3d> truth;
    std::mt19937 rng(3);
    std::uniform_int_distribution<int> hkl(-5, 5);

    // Sweep φ so enough (hkl,φ) pairs land on the Ewald sphere.
    for (int h = -4; h <= 4; ++h) {
        for (int k = -4; k <= 4; ++k) {
            for (int l = -4; l <= 4; ++l) {
                if (h == 0 && k == 0 && l == 0) {
                    continue;
                }
                const Eigen::Vector3d r0 =
                    amat * Eigen::Vector3d(h, k, l);
                for (int ip = 0; ip < 36; ++ip) {
                    const double phi = ip * 5.0;
                    double x_mm = 0.0;
                    double y_mm = 0.0;
                    if (!rlp_to_spot_mm(
                            r0, phi, geo.distance_mm, wave, geo.beam_x_mm,
                            geo.beam_y_mm, x_mm, y_mm
                        )) {
                        continue;
                    }
                    Spot s;
                    s.x_mm = x_mm;
                    s.y_mm = y_mm;
                    s.phi_start_deg = phi;
                    s.phi_end_deg = phi;
                    s.intensity = 100.0;
                    s.sigma = 1.0;
                    spots.spots.push_back(s);
                    truth.push_back(r0);
                    if (truth.size() >= 60u) {
                        goto enough;
                    }
                }
            }
        }
    }
enough:
    CHECK(truth.size() >= 20u);

    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> z;
    const std::size_t n =
        spots_to_rlp_xyz(spots, {}, geo, x, y, z);
    CHECK(n == truth.size());

    double max_err = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const Eigen::Vector3d got(
            static_cast<double>(x[i]),
            static_cast<double>(y[i]),
            static_cast<double>(z[i])
        );
        max_err = std::max(max_err, (got - truth[i]).norm());
    }
    std::cout << "rlp_roundtrip in-memory: n=" << n
              << " max_err=" << max_err << " Å⁻¹\n";
    CHECK_MSG(max_err < 1e-3, "in-memory RLP round-trip failed");

    // .spt write/read: still (start==end) → field 4 = phi_end = true φ.
    const auto tmp =
        std::filesystem::temp_directory_path() / "pl_rlp_roundtrip.spt";
    spots.write_spt(tmp);
    const SpotList loaded = SpotList::read_spt(tmp);
    std::error_code ec;
    std::filesystem::remove(tmp, ec);

    std::vector<float> xr;
    std::vector<float> yr;
    std::vector<float> zr;
    const std::size_t nr =
        spots_to_rlp_xyz(loaded, {}, geo, xr, yr, zr);
    CHECK(nr == truth.size());
    max_err = 0.0;
    for (std::size_t i = 0; i < nr; ++i) {
        const Eigen::Vector3d got(
            static_cast<double>(xr[i]),
            static_cast<double>(yr[i]),
            static_cast<double>(zr[i])
        );
        max_err = std::max(max_err, (got - truth[i]).norm());
    }
    std::cout << "rlp_roundtrip via .spt: n=" << nr
              << " max_err=" << max_err << " Å⁻¹\n";
    CHECK_MSG(max_err < 1e-3, ".spt RLP round-trip failed");

    (void)rng;
    return mosflm_test::finish("rlp_roundtrip");
}
