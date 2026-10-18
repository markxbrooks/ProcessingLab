#pragma once

#include "mosflm/diffraction/ProcessPipeline.hpp"
#include "mosflm/geometry/UnitCell.hpp"
#include "mosflm/image/DiffractionImage.hpp"

#include <Eigen/Dense>

#include <filesystem>
#include <iosfwd>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace mosflm::control {

/// Keyword state accumulated by a thin MOSFLM CONTROL subset (Phase D / D5).
struct ControlState {
    bool has_distance{false};
    bool has_wavelength{false};
    bool has_beam{false};
    bool has_pixel{false};
    bool has_cell{false};
    bool has_symmetry{false};
    bool has_matrix{false};
    bool has_mosaic{false};
    bool has_divergence{false};
    bool has_separation{false};

    double distance_mm{160.0};
    double wavelength_A{0.9790};
    double beam_x_mm{94.0};
    double beam_y_mm{94.7};
    double pixel_x_mm{0.0816};
    double pixel_y_mm{0.0816};
    double mosaic_deg{0.30};
    /// Beam divergence full-angles (degrees); MOSFLM ``DIVERGENCE`` / DIVH,DIVV.
    double divh_deg{0.02};
    double divv_deg{0.02};
    double separation_x_mm{0.57};
    double separation_y_mm{0.57};
    double dmin_A{1.45};

    UnitCell cell{78.19, 78.19, 78.19, 90.0, 90.0, 90.0};
    int symmetry_number{1};
    /// MOSFLM Orientation Matrix [A] (lab X,Y,Z rows), before C++ frame map.
    Eigen::Matrix3d amat{Eigen::Matrix3d::Identity()};

    std::filesystem::path directory{"."};
    std::string image_template;  ///< e.g. ``insulin_1_###.img``
    std::filesystem::path image_path;  ///< ``IMAGE`` keyword if set
    /// Output MTZ path (``HKLOUT`` keyword or CLI override).
    std::filesystem::path hklout{"process.mtz"};

    bool process_pending{false};
    int process_first{1};
    int process_last{1};
    double process_start_phi_deg{0.0};
    double process_angle_deg{1.0};

    bool exit_requested{false};
};

/// Outcome of one ``GO`` that ran PROCESS.
struct ControlGoResult {
    bool ran_process{false};
    diffraction::ProcessResult process{};
    crystal::CrystalMetadata crystal{};
    std::filesystem::path image_used{};
    std::string message;
};

/// Pure-C++ CONTROL driver for the known-indexed keyword subset.
///
/// Supported keywords: ``CELL``, ``SYMMETRY``, ``BEAM``, ``DISTANCE``,
/// ``WAVELENGTH``, ``PIXEL``, ``MATRIX`` / ``AMAT``, ``MOSAIC``,
/// ``DIVERGENCE``, ``SEPARATION``, ``DIRECTORY``, ``TEMPLATE``, ``IMAGE``,
/// ``HKLOUT``, ``PROCESS``, ``GO``, ``EXIT``. Unknown lines (``AUTOINDEX``,
/// ``GENFILE``, …) are ignored so a trimmed ``ipmosflm.in`` remains readable.
class ControlSession {
public:
    /// Parse one logical line (comments stripped). Throws on malformed
    /// known keywords.
    void parse_line(std::string_view line);

    /// Parse a whole script file.
    void parse_file(const std::filesystem::path& path);

    /// Execute ``GO``: if ``PROCESS`` is pending, load the image and run
    /// :class:`diffraction::ProcessPipeline`.
    ///
    /// ``image_override`` replaces ``IMAGE`` / ``TEMPLATE`` resolution (tests
    /// and CLI ``--image``).
    [[nodiscard]] ControlGoResult go(
        const std::optional<std::filesystem::path>& image_override = std::nullopt
    );

    /// Build a :struct:`diffraction::ProcessRequest` from keyword state + image.
    [[nodiscard]] diffraction::ProcessRequest make_process_request(
        const image::DiffractionImage& image
    ) const;

    /// Resolve the SMV path for the pending / last PROCESS frame.
    [[nodiscard]] std::filesystem::path resolve_image_path(
        const std::optional<std::filesystem::path>& image_override =
            std::nullopt
    ) const;

    [[nodiscard]] const ControlState& state() const { return state_; }
    [[nodiscard]] bool exit_requested() const { return state_.exit_requested; }

    /// Parse an input stream and run every ``GO`` until ``EXIT`` (or EOF).
    [[nodiscard]] std::vector<ControlGoResult> run_stream(
        std::istream& in,
        const std::optional<std::filesystem::path>& image_override =
            std::nullopt
    );

    /// Parse file and run every ``GO`` until ``EXIT`` (or EOF).
    [[nodiscard]] std::vector<ControlGoResult> run_script(
        const std::filesystem::path& script,
        const std::optional<std::filesystem::path>& image_override =
            std::nullopt
    );

private:
    ControlState state_;
    /// Numbers still needed to finish a multi-line ``MATRIX``.
    int matrix_needed_{0};
    std::vector<double> matrix_buf_;

    void apply_keyword(const std::string& key, const std::vector<std::string>& args);
    void ingest_matrix_numbers(const std::vector<double>& vals);
    void finish_matrix_if_ready();
};

}  // namespace mosflm::control
