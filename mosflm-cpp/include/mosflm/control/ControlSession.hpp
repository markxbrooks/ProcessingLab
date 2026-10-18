#pragma once

#include "mosflm/diffraction/ProcessPipeline.hpp"
#include "mosflm/geometry/UnitCell.hpp"
#include "mosflm/image/DiffractionImage.hpp"

#include <Eigen/Dense>

#include <filesystem>
#include <functional>
#include <iosfwd>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <vector>

namespace mosflm::control {

/// AUTOINDEX settings accumulated from CONTROL (Phase E / E1).
struct AutoIndexControlSettings {
    std::filesystem::path spot_file;
    float min_isigi{20.0f};
    float hkldev{0.3f};
    int nvect{30};
    int nvref{30};
    int prefer_character{0};
    std::string prefer_system;
    int image_number{1};
};

/// Keyword state accumulated by MOSFLM CONTROL (Phase D / E).
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
    /// Default 0 unless ``DIVERGENCE`` is set (matches known-indexed insulin).
    double divh_deg{0.0};
    double divv_deg{0.0};
    double separation_x_mm{0.57};
    double separation_y_mm{0.57};
    double dmin_A{1.45};
    double dmax_A{200.0};
    /// Predictor φ sample step (degrees); 0.025 matches known_indexed_insulin.
    double phi_step_deg{0.025};
    /// Oracle-refined detector fiducials (insulin PROCESS / geometry.txt).
    double yscale{0.9921};
    double tilt_mrad{-83.0};
    double twist_mrad{-140.0};
    bool use_rotmat2{true};

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

    /// Standard-profile detector tiling (MOSFLM NMASKS / XLINES×YLINES).
    int nmasks_x{5};
    int nmasks_y{5};
    /// Opt-in Stage 2 post-refinement on PROCESS.
    bool enable_postrefinement{false};
    /// Multi-image NADD window (1 = single image; >1 reserved for joint refine).
    int nadd{1};

    bool autoindex_pending{false};
    AutoIndexControlSettings autoindex{};
    /// Filled by AUTOINDEX runner for DiffractionLab / machine-readable dump.
    std::string autoindex_system;
    int autoindex_character{0};
    std::size_t autoindex_n_spots{0};

    /// Multi-lattice: additional orientation matrices (first is ``amat``).
    std::vector<Eigen::Matrix3d> extra_lattices;

    bool strategy_pending{false};
    bool strategy_done{false};
    double strategy_osc_range_deg{90.0};
    double strategy_completeness_est{0.0};

    bool done_run{false};
    bool exit_requested{false};
};

/// Outcome of one ``GO`` / ``RUN``.
struct ControlGoResult {
    bool ran_process{false};
    bool ran_autoindex{false};
    bool ran_strategy{false};
    diffraction::ProcessResult process{};
    crystal::CrystalMetadata crystal{};
    std::filesystem::path image_used{};
    std::string message;
};

/// Hook that runs AUTOINDEX and writes CELL/MATRIX into ``state``.
using AutoIndexRunner = std::function<void(ControlState& state)>;

/// Pure-C++ CONTROL driver (Phase D subset + Phase E extensions).
///
/// Supported: geometry keywords, ``PROCESS``, ``GO``/``RUN``, ``EXIT``,
/// ``AUTOINDEX`` (when a runner is registered), ``STRATEGY``, ``NMASKS``,
/// ``POSTREFINE``, ``NADD``, ``MULTILATTICE``. Unimplemented major verbs
/// (``REFIX``, …) throw. Other unknown keywords warn once then continue.
class ControlSession {
public:
    /// Register DPS AUTOINDEX implementation (from ``mosflm_indexing``).
    static void register_autoindex_runner(AutoIndexRunner runner);

    /// Parse one logical line (comments stripped). Throws on malformed
    /// known keywords or unimplemented major verbs.
    void parse_line(std::string_view line);

    /// Parse a whole script file.
    void parse_file(const std::filesystem::path& path);

    /// Execute ``GO``/``RUN``: AUTOINDEX, STRATEGY, and/or PROCESS.
    [[nodiscard]] ControlGoResult go(
        const std::optional<std::filesystem::path>& image_override = std::nullopt
    );

    [[nodiscard]] diffraction::ProcessRequest make_process_request(
        const image::DiffractionImage& image
    ) const;

    [[nodiscard]] std::filesystem::path resolve_image_path(
        const std::optional<std::filesystem::path>& image_override =
            std::nullopt
    ) const;

    [[nodiscard]] const ControlState& state() const { return state_; }
    [[nodiscard]] bool exit_requested() const { return state_.exit_requested; }
    [[nodiscard]] bool done_run() const { return state_.done_run; }
    [[nodiscard]] bool strategy_done() const { return state_.strategy_done; }

    /// True when EOF should emit Fortran-style no-RUN warning.
    [[nodiscard]] bool needs_no_run_warning() const {
        return !state_.done_run && !state_.strategy_done;
    }

    [[nodiscard]] std::vector<ControlGoResult> run_stream(
        std::istream& in,
        const std::optional<std::filesystem::path>& image_override =
            std::nullopt
    );

    [[nodiscard]] std::vector<ControlGoResult> run_script(
        const std::filesystem::path& script,
        const std::optional<std::filesystem::path>& image_override =
            std::nullopt
    );

private:
    ControlState state_;
    int matrix_needed_{0};
    std::vector<double> matrix_buf_;
    std::set<std::string> warned_unknown_;

    void apply_keyword(const std::string& key, const std::vector<std::string>& args);
    void ingest_matrix_numbers(const std::vector<double>& vals);
    void finish_matrix_if_ready();
    void parse_autoindex_args(const std::vector<std::string>& args);
};

}  // namespace mosflm::control
