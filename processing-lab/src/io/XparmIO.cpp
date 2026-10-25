#include "mosflm/io/XparmIO.hpp"

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/geometry/LabFrame.hpp"
#include "mosflm/io/ImportedSweepGeometry.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mosflm::io {
namespace {

std::vector<std::string> split_ws(const std::string& line) {
    std::istringstream iss(line);
    std::vector<std::string> out;
    std::string tok;
    while (iss >> tok) {
        out.push_back(tok);
    }
    return out;
}

bool nearly_orthonormal_identity(const Eigen::Matrix3d& A) {
    const double nrm = A.norm();
    if (nrm < 1.2 || nrm > 2.2) {
        return false;
    }
    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    return (A - I).norm() < 0.15 || (A + I).norm() < 0.15 ||
           (A.transpose() * A - I).norm() < 0.2;
}

/// Map from the XDS lab frame to ProcessingLab's lab frame. Rows are the
/// new x, y, z: z is the beam, x the rotation axis made normal to the beam.
/// ProcessingLab predicts with the beam along +z and rotation about x, so
/// XDS vectors must be brought into that frame; otherwise a tilted beam or
/// axis (0.5–0.6° on many beamlines) shifts φ by ±0.5° and positions.
///
/// x is negated: with the detector omega of an axis along +fast (CCOM 90°,
/// as for Pilatus "X, CW" and the thaumatin Eiger) ProcessingLab's detector
/// mapping has the opposite x handedness to XDS, and predictions come out
/// mirrored about the beam (XDET = 2·beam − XCAL). The mirror commutes with
/// rotation about x and leaves the beam unchanged, so φ and the Miller
/// indices are unaffected (checked against XDS INTEGRATE.HKL: 1.6 px, 0.02°).
Eigen::Matrix3d xds_to_lab_frame(const double beam[3], const double axis[3]) {
    const Eigen::Vector3d z = Eigen::Vector3d(beam[0], beam[1], beam[2]).normalized();
    Eigen::Vector3d x(axis[0], axis[1], axis[2]);
    x = (x - x.dot(z) * z).normalized();
    const Eigen::Vector3d y = z.cross(x);
    Eigen::Matrix3d R;
    R.row(0) = -x;
    R.row(1) = y;
    R.row(2) = z;
    return R;
}

/// Build Mosflm AMAT from XDS real-space unit-cell axes (lab Å).
Eigen::Matrix3d mosflm_amat_from_xds_real_axes(
    const Eigen::Vector3d& ax,
    const Eigen::Vector3d& bx,
    const Eigen::Vector3d& cx,
    double wavelength_A
) {
    Eigen::Matrix3d real;
    real.col(0) = ax;
    real.col(1) = bx;
    real.col(2) = cx;
    const double det = real.determinant();
    if (!(std::fabs(det) > 1.0)) {
        return Eigen::Matrix3d::Zero();
    }
    // Reciprocal basis (Å⁻¹) in the same lab frame as the real axes.
    const Eigen::Matrix3d recip = real.inverse().transpose();
    // Treat XDS/Xia2 lab as the C++ lab; convert to MOSFLM AMAT.
    return cpp_ub_to_mosflm_amat(recip, wavelength_A);
}

bool try_parse_xds_classic(
    const std::vector<std::vector<std::string>>& rows,
    XparmGeometry& g
) {
    // Locate space-group + cell row (7 tokens).
    std::size_t icell = rows.size();
    for (std::size_t i = 0; i < rows.size(); ++i) {
        if (rows[i].size() < 7) {
            continue;
        }
        try {
            const int sg = std::stoi(rows[i][0]);
            const double a = std::stod(rows[i][1]);
            const double b = std::stod(rows[i][2]);
            const double c = std::stod(rows[i][3]);
            const double al = std::stod(rows[i][4]);
            const double be = std::stod(rows[i][5]);
            const double ga = std::stod(rows[i][6]);
            if (a > 1.0 && b > 1.0 && c > 1.0 && al > 10.0 && be > 10.0 &&
                ga > 10.0 && al < 180.0 && sg >= 1 && sg <= 230) {
                g.spacegroup_number = sg;
                g.cell = UnitCell{a, b, c, al, be, ga};
                icell = i;
                break;
            }
        } catch (...) {
        }
    }
    if (icell >= rows.size()) {
        return false;
    }

    // Wavelength + beam: exactly 4 tokens, wave in [0.3,3]. (The 6-token
    // "starting frame, angle, osc, axis" row before it starts with the frame
    // number, e.g. 1, which would otherwise be read as λ = 1 Å.)
    for (const auto& r : rows) {
        if (r.size() != 4) {
            continue;
        }
        try {
            const double wave = std::stod(r[0]);
            if (wave > 0.3 && wave < 3.0) {
                g.wavelength_A = wave;
                g.beam_dir[0] = std::stod(r[1]);
                g.beam_dir[1] = std::stod(r[2]);
                g.beam_dir[2] = std::stod(r[3]);
                break;
            }
        } catch (...) {
        }
    }

    // Starting frame / osc / rotation axis (6 tokens).
    for (const auto& r : rows) {
        if (r.size() < 6) {
            continue;
        }
        try {
            const double sf = std::stod(r[0]);
            const double sa = std::stod(r[1]);
            const double osc = std::stod(r[2]);
            if (sf >= 1.0 && sf < 100000.0 && osc > 0.001 && osc < 5.0) {
                g.starting_frame = sf;
                g.starting_angle_deg = sa;
                g.oscillation_range_deg = osc;
                g.rot_axis[0] = std::stod(r[3]);
                g.rot_axis[1] = std::stod(r[4]);
                g.rot_axis[2] = std::stod(r[5]);
                break;
            }
        } catch (...) {
        }
    }

    // Three real-space cell-axis rows immediately after the cell line.
    if (icell + 3 >= rows.size()) {
        return false;
    }
    Eigen::Vector3d ax, bx, cx;
    bool axes_ok = true;
    try {
        for (int k = 0; k < 3; ++k) {
            const auto& r = rows[icell + 1 + static_cast<std::size_t>(k)];
            if (r.size() != 3) {
                axes_ok = false;
                break;
            }
            Eigen::Vector3d v(
                std::stod(r[0]), std::stod(r[1]), std::stod(r[2])
            );
            const double len = v.norm();
            // Real-space edges are tens–hundreds of Å (not unit detector axes).
            if (!(len > 5.0 && len < 2000.0)) {
                axes_ok = false;
                break;
            }
            if (k == 0) {
                ax = v;
            } else if (k == 1) {
                bx = v;
            } else {
                cx = v;
            }
        }
    } catch (...) {
        axes_ok = false;
    }
    if (axes_ok) {
        const Eigen::Matrix3d R = xds_to_lab_frame(g.beam_dir, g.rot_axis);
        const Eigen::Matrix3d A = mosflm_amat_from_xds_real_axes(
            R * ax, R * bx, R * cx, g.wavelength_A > 0.1 ? g.wavelength_A : 1.0
        );
        if (A.norm() > 1e-6) {
            g.amat = A;
            g.has_amat = true;
        }
    }

    // Detector: "segment NX NY QX QY" then "ORGX ORGY DISTANCE".
    for (std::size_t i = 0; i + 1 < rows.size(); ++i) {
        const auto& a = rows[i];
        const auto& b = rows[i + 1];
        if (a.size() == 5 && b.size() >= 3) {
            try {
                const int nx = std::stoi(a[1]);
                const int ny = std::stoi(a[2]);
                const double qx = std::stod(a[3]);
                const double qy = std::stod(a[4]);
                const double ox = std::stod(b[0]);
                const double oy = std::stod(b[1]);
                const double dist = std::stod(b[2]);
                if (nx > 10 && ny > 10 && qx > 0.01 && qx < 1.0 &&
                    qy > 0.01 && qy < 1.0 && dist > 10.0 && ox > 10.0 &&
                    oy > 10.0) {
                    g.nfast = nx;
                    g.nslow = ny;
                    g.qx_mm = qx;
                    g.qy_mm = qy;
                    g.orgx_px = ox;
                    g.orgy_px = oy;
                    g.distance_mm = dist;
                    // Detector X and Y axis rows follow ORGX ORGY F.
                    for (int k = 0; k < 2 && i + 2 + k < rows.size(); ++k) {
                        const auto& d = rows[i + 2 + static_cast<std::size_t>(k)];
                        if (d.size() == 3) {
                            double* dst = (k == 0) ? g.det_x : g.det_y;
                            for (int c = 0; c < 3; ++c) {
                                dst[c] = std::stod(d[static_cast<std::size_t>(c)]);
                            }
                        }
                    }
                    break;
                }
            } catch (...) {
            }
        }
    }

    // Direct-beam pixel: intersect the incident beam with the detector plane
    // (normal N = X × Y at distance F from the crystal).
    {
        const Eigen::Vector3d s0 =
            Eigen::Vector3d(g.beam_dir[0], g.beam_dir[1], g.beam_dir[2]).normalized();
        const Eigen::Vector3d dx(g.det_x[0], g.det_x[1], g.det_x[2]);
        const Eigen::Vector3d dy(g.det_y[0], g.det_y[1], g.det_y[2]);
        const Eigen::Vector3d n = dx.cross(dy).normalized();
        const double cosn = s0.dot(n);
        if (cosn > 0.5 && g.qx_mm > 0.0 && g.qy_mm > 0.0 && g.orgx_px > 0.0) {
            const Eigen::Vector3d p = s0 * (g.distance_mm / cosn);
            g.beam_fast_px = g.orgx_px + p.dot(dx) / g.qx_mm;
            g.beam_slow_px = g.orgy_px + p.dot(dy) / g.qy_mm;
            g.has_direct_beam = true;
            // The detector normal in the frame of xds_to_lab_frame (rows −x,
            // y = z × x, z) gives TILT / TWIST: n = (−TWIST, −TILT, 1000)/1000
            // (checked on 76c1b, where refining them from zero converged on
            // the same values).
            const Eigen::Vector3d nl = xds_to_lab_frame(g.beam_dir, g.rot_axis) * n;
            g.tilt_mrad = -1000.0 * nl.y() / nl.z();
            g.twist_mrad = -1000.0 * nl.x() / nl.z();
            g.has_tilt_twist = std::fabs(g.tilt_mrad) > 1e-3 || std::fabs(g.twist_mrad) > 1e-3;
        }
    }

    return g.cell.a > 1.0 && g.has_amat;
}

}  // namespace

XparmGeometry read_xparm(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("read_xparm: cannot open " + path.string());
    }
    XparmGeometry g;
    std::string line;
    std::vector<std::vector<std::string>> rows;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '!' || line[0] == '#') {
            continue;
        }
        // Skip the "XPARM.XDS VERSION ..." banner.
        if (line.find("XPARM") != std::string::npos &&
            line.find("VERSION") != std::string::npos) {
            continue;
        }
        auto toks = split_ws(line);
        if (!toks.empty()) {
            rows.push_back(std::move(toks));
        }
    }

    // Prefer a full Xia2/XDS classic layout (cell axes → Mosflm AMAT).
    if (try_parse_xds_classic(rows, g)) {
        return g;
    }

    // Fallback: ProcessingLab bridge GXPARM / loose keyword scan.
    for (const auto& r : rows) {
        if (r.size() < 7) {
            continue;
        }
        try {
            const int sg = std::stoi(r[0]);
            const double a = std::stod(r[1]);
            const double b = std::stod(r[2]);
            const double c = std::stod(r[3]);
            const double al = std::stod(r[4]);
            const double be = std::stod(r[5]);
            const double ga = std::stod(r[6]);
            if (a > 1.0 && b > 1.0 && c > 1.0 && al > 10.0 && be > 10.0 &&
                ga > 10.0 && al < 180.0 && sg >= 1 && sg <= 230) {
                g.spacegroup_number = sg;
                g.cell = UnitCell{a, b, c, al, be, ga};
                break;
            }
        } catch (...) {
        }
    }

    for (const auto& r : rows) {
        if (r.size() < 7) {
            continue;
        }
        try {
            const int nx = std::stoi(r[0]);
            const int ny = std::stoi(r[1]);
            const double ox = std::stod(r[2]);
            const double oy = std::stod(r[3]);
            const double dist = std::stod(r[4]);
            const double qx = std::stod(r[5]);
            const double qy = std::stod(r[6]);
            // ProcessingLab's bridge GXPARM writes NX NY as 0 0 (unknown).
            const bool size_ok = (nx == 0 && ny == 0) || (nx > 10 && ny > 10);
            if (size_ok && dist > 10.0 && qx > 0.01 && qx < 1.0 &&
                qy > 0.01 && qy < 1.0) {
                g.nfast = nx;
                g.nslow = ny;
                g.orgx_px = ox;
                g.orgy_px = oy;
                g.distance_mm = dist;
                g.qx_mm = qx;
                g.qy_mm = qy;
            }
        } catch (...) {
        }
    }

    if (g.distance_mm <= 10.0) {
        for (const auto& r : rows) {
            if (r.size() < 5) {
                continue;
            }
            try {
                const double ox = std::stod(r[0]);
                const double oy = std::stod(r[1]);
                const double dist = std::stod(r[2]);
                const double qx = std::stod(r[3]);
                const double qy = std::stod(r[4]);
                if (dist > 10.0 && qx > 0.01 && qx < 1.0 && qy > 0.01 &&
                    qy < 1.0 && ox > 10.0 && oy > 10.0) {
                    g.orgx_px = ox;
                    g.orgy_px = oy;
                    g.distance_mm = dist;
                    g.qx_mm = qx;
                    g.qy_mm = qy;
                }
            } catch (...) {
            }
        }
    }

    // Mosflm-scale AMAT only (reject detector-axis identity blocks).
    for (std::size_t i = 0; i + 2 < rows.size(); ++i) {
        if (rows[i].size() != 3 || rows[i + 1].size() != 3 ||
            rows[i + 2].size() != 3) {
            continue;
        }
        try {
            Eigen::Matrix3d A;
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) {
                    A(r, c) = std::stod(rows[i + static_cast<std::size_t>(r)]
                                            [static_cast<std::size_t>(c)]);
                }
            }
            const double nrm = A.norm();
            if (nrm > 1e-6 && nrm < 10.0 && !nearly_orthonormal_identity(A)) {
                g.amat = A;
                g.has_amat = true;
                break;
            }
        } catch (...) {
        }
    }

    for (const auto& r : rows) {
        if (r.size() < 4) {
            continue;
        }
        try {
            const double wave = std::stod(r[0]);
            if (wave > 0.3 && wave < 3.0) {
                g.wavelength_A = wave;
                g.beam_dir[0] = std::stod(r[1]);
                g.beam_dir[1] = std::stod(r[2]);
                g.beam_dir[2] = std::stod(r[3]);
            }
        } catch (...) {
        }
    }

    return g;
}

void write_xparm(
    const std::filesystem::path& path,
    const XparmGeometry& g
) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("write_xparm: cannot open " + path.string());
    }
    out << std::setprecision(10);
    out << "! GXPARM.XDS written by ProcessingLab (bridge format)\n";
    out << g.starting_frame << ' ' << g.starting_angle_deg << ' '
        << g.oscillation_range_deg << ' ' << g.rot_axis[0] << ' '
        << g.rot_axis[1] << ' ' << g.rot_axis[2] << '\n';
    out << g.wavelength_A << ' ' << g.beam_dir[0] << ' ' << g.beam_dir[1]
        << ' ' << g.beam_dir[2] << '\n';
    out << g.spacegroup_number << ' ' << g.cell.a << ' ' << g.cell.b << ' '
        << g.cell.c << ' ' << g.cell.alpha_deg << ' ' << g.cell.beta_deg
        << ' ' << g.cell.gamma_deg << '\n';
    if (g.has_amat) {
        for (int r = 0; r < 3; ++r) {
            out << g.amat(r, 0) << ' ' << g.amat(r, 1) << ' ' << g.amat(r, 2)
                << '\n';
        }
    } else {
        out << "0.01 0 0\n0 0.01 0\n0 0 0.01\n";
    }
    out << g.nfast << ' ' << g.nslow << ' ' << g.orgx_px << ' ' << g.orgy_px
        << ' ' << g.distance_mm << ' ' << g.qx_mm << ' ' << g.qy_mm << '\n';
    out << g.det_x[0] << ' ' << g.det_x[1] << ' ' << g.det_x[2] << '\n';
    out << g.det_y[0] << ' ' << g.det_y[1] << ' ' << g.det_y[2] << '\n';
}

namespace {

/// Where the beam meets the detector (mm), if the XPARM gives it.
std::optional<std::pair<double, double>> xparm_beam_mm(const XparmGeometry& g) {
    if (!(g.qx_mm > 0.0 && g.qy_mm > 0.0)) {
        return std::nullopt;
    }
    if (g.has_direct_beam) {
        // Classic XDS: BEAM is where the beam meets the detector (fast,
        // slow), not ORGX/ORGY (the foot of the detector normal).
        return std::pair{g.beam_fast_px * g.qx_mm, g.beam_slow_px * g.qy_mm};
    }
    if (g.orgx_px > 0.0 && g.orgy_px > 0.0) {
        // ProcessingLab bridge GXPARM stores its own BEAM as ORGX/ORGY.
        return std::pair{g.orgx_px * g.qx_mm, g.orgy_px * g.qy_mm};
    }
    return std::nullopt;
}

}  // namespace

ImportedSweepGeometry imported_from_xparm(const XparmGeometry& g) {
    ImportedSweepGeometry im;
    if (g.cell.a > 1.0) {
        im.has_cell = true;
        im.cell_a = g.cell.a;
        im.cell_b = g.cell.b;
        im.cell_c = g.cell.c;
        im.cell_alpha = g.cell.alpha_deg;
        im.cell_beta = g.cell.beta_deg;
        im.cell_gamma = g.cell.gamma_deg;
    }
    if (g.spacegroup_number >= 1 && g.spacegroup_number <= 230) {
        im.has_symmetry = true;
        im.symmetry_number = g.spacegroup_number;
    }
    if (g.wavelength_A > 0.1) {
        im.has_wavelength = true;
        im.wavelength_A = g.wavelength_A;
    }
    if (g.distance_mm > 1.0) {
        im.has_distance = true;
        im.distance_mm = g.distance_mm;
    }
    if (g.qx_mm > 0.0 && g.qy_mm > 0.0) {
        im.has_pixel = true;
        im.pixel_x_mm = g.qx_mm;
        im.pixel_y_mm = g.qy_mm;
    }
    if (const auto b = xparm_beam_mm(g)) {
        im.has_beam = true;
        im.beam_x_mm = b->first;
        im.beam_y_mm = b->second;
    }
    return im;
}

std::vector<std::string> xparm_orientation_lines(const XparmGeometry& g) {
    if (!g.has_amat) {
        return {};
    }
    std::ostringstream m;
    m << std::setprecision(12) << "MATRIX";
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            m << ' ' << g.amat(r, c);
        }
    }
    return {
        // The AMAT is built in the XDS lab frame (beam +z, axis +x), whose
        // handedness matches ProcessingLab .spt indexing, not MOSFLM's.
        "INVERTX OFF",
        // Keep the XDS-refined detector and cell; refine the orientation
        // only, as XDS INTEGRATE does. Distance and cell together are
        // degenerate (both scale the pattern; 385 → 390 mm on 76c1b), and
        // per-period cell refinement wandered b 185–189 Å late in that sweep.
        "REFINE ORIENTATION",
        m.str(),
    };
}

void apply_xparm(
    control::ControlSession& session,
    const XparmGeometry& g
) {
    std::ostringstream line;
    line << std::setprecision(12);
    const auto emit = [&] {
        session.parse_line(line.str());
        line.str("");
        line.clear();
    };
    if (g.cell.a > 1.0) {
        line << "CELL " << g.cell.a << ' ' << g.cell.b << ' ' << g.cell.c << ' '
             << g.cell.alpha_deg << ' ' << g.cell.beta_deg << ' '
             << g.cell.gamma_deg;
        emit();
    }
    if (g.spacegroup_number >= 1 && g.spacegroup_number <= 230) {
        session.parse_line(
            "SYMMETRY " + std::to_string(g.spacegroup_number)
        );
    }
    if (g.wavelength_A > 0.1) {
        line << "WAVELENGTH " << g.wavelength_A;
        emit();
    }
    if (g.distance_mm > 1.0) {
        line << "DISTANCE " << g.distance_mm;
        emit();
    }
    if (g.qx_mm > 0.0 && g.qy_mm > 0.0) {
        line << "PIXEL " << g.qx_mm << ' ' << g.qy_mm;
        emit();
    }
    if (const auto b = xparm_beam_mm(g)) {
        line << "BEAM " << b->first << ' ' << b->second;
        emit();
    }
    if (g.has_tilt_twist) {
        line << "TILT " << g.tilt_mrad;
        emit();
        line << "TWIST " << g.twist_mrad;
        emit();
    }
    for (const auto& l : xparm_orientation_lines(g)) {
        session.parse_line(l);
    }
}

XparmGeometry xparm_from_control(const control::ControlSession& session) {
    const auto& st = session.state();
    XparmGeometry g;
    g.cell = st.cell;
    g.spacegroup_number = st.has_symmetry ? st.symmetry_number : 1;
    g.wavelength_A = st.wavelength_A;
    g.distance_mm = st.distance_mm;
    g.qx_mm = st.pixel_x_mm;
    g.qy_mm = st.pixel_y_mm;
    if (st.pixel_x_mm > 0.0 && st.pixel_y_mm > 0.0) {
        g.orgx_px = st.beam_x_mm / st.pixel_x_mm;
        g.orgy_px = st.beam_y_mm / st.pixel_y_mm;
    }
    if (st.has_matrix) {
        g.amat = st.amat;
        g.has_amat = true;
    }
    return g;
}

}  // namespace mosflm::io
