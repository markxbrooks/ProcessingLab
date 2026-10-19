#include "mosflm/io/XparmIO.hpp"

#include "mosflm/control/ControlSession.hpp"
#include "mosflm/geometry/LabFrame.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
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

    // Wavelength + beam (4 tokens, wave in [0.3,3]).
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
        const Eigen::Matrix3d A = mosflm_amat_from_xds_real_axes(
            ax, bx, cx, g.wavelength_A > 0.1 ? g.wavelength_A : 1.0
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
                    break;
                }
            } catch (...) {
            }
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
            if (nx > 10 && ny > 10 && dist > 10.0 && qx > 0.01 && qx < 1.0 &&
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

void apply_xparm(
    control::ControlSession& session,
    const XparmGeometry& g
) {
    std::ostringstream line;
    line << std::setprecision(12);
    if (g.cell.a > 1.0) {
        line.str("");
        line.clear();
        line << "CELL " << g.cell.a << ' ' << g.cell.b << ' ' << g.cell.c << ' '
             << g.cell.alpha_deg << ' ' << g.cell.beta_deg << ' '
             << g.cell.gamma_deg;
        session.parse_line(line.str());
    }
    if (g.spacegroup_number >= 1 && g.spacegroup_number <= 230) {
        session.parse_line(
            "SYMMETRY " + std::to_string(g.spacegroup_number)
        );
    }
    if (g.wavelength_A > 0.1) {
        line.str("");
        line.clear();
        line << "WAVELENGTH " << g.wavelength_A;
        session.parse_line(line.str());
    }
    if (g.distance_mm > 1.0) {
        line.str("");
        line.clear();
        line << "DISTANCE " << g.distance_mm;
        session.parse_line(line.str());
    }
    if (g.qx_mm > 0.0 && g.qy_mm > 0.0) {
        line.str("");
        line.clear();
        line << "PIXEL " << g.qx_mm << ' ' << g.qy_mm;
        session.parse_line(line.str());
        if (g.orgx_px > 0.0 && g.orgy_px > 0.0) {
            line.str("");
            line.clear();
            line << "BEAM " << (g.orgx_px * g.qx_mm) << ' '
                 << (g.orgy_px * g.qy_mm);
            session.parse_line(line.str());
        }
    }
    if (g.has_amat) {
        line.str("");
        line.clear();
        line << "MATRIX";
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                line << ' ' << g.amat(r, c);
            }
        }
        session.parse_line(line.str());
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
