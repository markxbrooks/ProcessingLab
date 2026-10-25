#pragma once

#include "mosflm/geometry/UnitCell.hpp"

#include <Eigen/Dense>

#include <filesystem>
#include <string>
#include <vector>

namespace mosflm::control {
class ControlSession;
}

namespace mosflm::io {

struct ImportedSweepGeometry;

/// Minimal XDS ``XPARM.XDS`` / ``GXPARM.XDS`` geometry snapshot.
struct XparmGeometry {
    int spacegroup_number{1};
    UnitCell cell{};
    /// Incident beam direction (XDS lab, unit vector).
    double beam_dir[3]{0.0, 0.0, 1.0};
    double wavelength_A{1.0};
    /// Rotation axis (unit vector).
    double rot_axis[3]{1.0, 0.0, 0.0};
    /// Starting oscillation angle (deg) and range per image.
    double starting_angle_deg{0.0};
    double oscillation_range_deg{0.1};
    double starting_frame{1.0};
    /// Detector: orgx/orgy (px), distance (mm), QX/QY (mm/px).
    double orgx_px{0.0};
    double orgy_px{0.0};
    double distance_mm{100.0};
    double qx_mm{0.075};
    double qy_mm{0.075};
    int nfast{0};
    int nslow{0};
    /// Detector X/Y axis directions (XDS lab). Identity defaults.
    double det_x[3]{1.0, 0.0, 0.0};
    double det_y[3]{0.0, 1.0, 0.0};
    /// Where the incident beam meets the detector (px, fast / slow), from a
    /// classic XDS file. ORGX/ORGY is the foot of the detector normal, which
    /// differs from the beam position when the beam is not normal to the
    /// detector (e.g. 22 px on a 0.56° beam tilt at 385 mm).
    bool has_direct_beam{false};
    double beam_fast_px{0.0};
    double beam_slow_px{0.0};
    /// Detector TILT / TWIST (mrad) in ProcessingLab's frame (beam +z, axis
    /// +x), from the detector normal there. Non-zero
    /// when XDS's beam is not along the detector normal (76c1b: a 10 mrad
    /// beam tilt left predictions 2 px out at the detector edge).
    bool has_tilt_twist{false};
    double tilt_mrad{0.0};
    double twist_mrad{0.0};
    /// Orientation: 3×3 A matrix (MOSFLM-like rows) if known. From a classic
    /// XDS file it is expressed in the frame with the beam along +z and the
    /// rotation axis along +x (ProcessingLab's lab frame).
    bool has_amat{false};
    Eigen::Matrix3d amat{Eigen::Matrix3d::Identity()};
};

/// Read a text XPARM.XDS / GXPARM.XDS (new-style multi-line layout).
[[nodiscard]] XparmGeometry read_xparm(const std::filesystem::path& path);

/// Write GXPARM.XDS-style text (compatible enough for Xia2/XDS bridges).
void write_xparm(
    const std::filesystem::path& path,
    const XparmGeometry& g
);

/// Cell, symmetry, wavelength, distance, pixel and beam from an XPARM, as
/// import fields (for per-field provenance; orientation not included).
[[nodiscard]] ImportedSweepGeometry imported_from_xparm(const XparmGeometry& g);

/// CONTROL lines that set the XPARM orientation (``INVERTX OFF``,
/// ``REFINE ORIENTATION``, ``MATRIX …``); empty if it has none.
[[nodiscard]] std::vector<std::string> xparm_orientation_lines(const XparmGeometry& g);

/// Apply XPARM cell / detector / matrix into CONTROL (does not set IMAGE).
void apply_xparm(
    control::ControlSession& session,
    const XparmGeometry& g
);

/// Snapshot current CONTROL geometry into an XPARM struct.
[[nodiscard]] XparmGeometry xparm_from_control(
    const control::ControlSession& session
);

}  // namespace mosflm::io
