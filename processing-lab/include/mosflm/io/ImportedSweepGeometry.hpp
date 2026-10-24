#pragma once

#include "mosflm/control/ControlSession.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace mosflm::io {

/// Geometry / sweep fields extracted from Xia2 ``.xinfo`` or ``XDS.INP``.
///
/// Optional fields use ``has_*`` flags. :meth:`to_control_lines` emits MOSFLM
/// CONTROL keywords for :class:`control::ControlSession`.
struct ImportedSweepGeometry {
    bool has_wavelength{false};
    bool has_distance{false};
    bool has_pixel{false};
    bool has_beam{false};
    bool has_cell{false};
    bool has_symmetry{false};
    bool has_process_range{false};
    bool has_oscillation{false};
    bool has_directory{false};
    bool has_image{false};
    bool has_template{false};

    double wavelength_A{0.0};
    double distance_mm{0.0};
    double pixel_x_mm{0.0};
    double pixel_y_mm{0.0};
    double beam_x_mm{0.0};
    double beam_y_mm{0.0};
    double cell_a{0.0};
    double cell_b{0.0};
    double cell_c{0.0};
    double cell_alpha{90.0};
    double cell_beta{90.0};
    double cell_gamma{90.0};
    int symmetry_number{1};
    int process_first{1};
    int process_last{1};
    double oscillation_deg{0.0};
    std::filesystem::path directory;
    std::filesystem::path image_path;
    std::string image_template;
    /// Untrusted regions and excluded resolution as CONTROL lines
    /// (``UNTRUSTED …``, ``EXCLUDE …``, ``DEFPIX VALUE_RANGE …``), from
    /// XDS.INP ``UNTRUSTED_*=`` etc. or ``run --untrusted``.
    std::vector<std::string> trusted_region_control;

    /// Overlay non-empty fields from ``other`` onto ``*this`` (other wins;
    /// ``trusted_region_control`` is appended).
    void merge_from(const ImportedSweepGeometry& other);

    /// Emit CONTROL lines for fields that are set.
    [[nodiscard]] std::vector<std::string> to_control_lines() const;
};

/// Apply geometry by feeding :meth:`ImportedSweepGeometry::to_control_lines`
/// into ``session.parse_line``.
void apply_imported_geometry(
    control::ControlSession& session,
    const ImportedSweepGeometry& geometry
);

/// Replace the last run of digits before the file extension with ``#`` of the
/// same width (e.g. ``sweep_0001.cbf`` → ``sweep_####.cbf``).
[[nodiscard]] std::string template_from_image_name(const std::string& filename);

/// Convert XDS ``?`` placeholders in a template basename to MOSFLM ``#``.
[[nodiscard]] std::string xds_question_marks_to_hash(const std::string& name);

}  // namespace mosflm::io
