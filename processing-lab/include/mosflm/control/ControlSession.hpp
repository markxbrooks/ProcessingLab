#pragma once

#include "mosflm/diffraction/ProcessPipeline.hpp"
#include "mosflm/geometry/UnitCell.hpp"
#include "mosflm/image/DiffractionImage.hpp"
#include "mosflm/util/OpenMP.hpp"

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
    /// XDS IDXREF knobs (Δ-vector path).
    float sepmin{6.0f};
    float cluster_radius{2.0f};
    bool use_diff_vector{true};
    /// Native DIALS-style 3D-FFT path (CONTROL ``FFT`` / ``NOFFT``).
    bool use_fft3d{true};
    /// Cambridge DPS path (off with ``FFTONLY``).
    bool use_dps{true};
    /// Geometry LS after INDEXFIT (Kabsch §2.8-lite). Off with ``NOIDXREF``.
    bool use_idxref_refine{true};
    float min_indexed_frac{0.70f};
};

/// FINDSPOTS settings (DiffractionLab-style peak finder → ``.spt``).
struct FindSpotsControlSettings {
    std::filesystem::path spot_file{"spots.spt"};
    int image_number{1};
    /// Optional multi-image COLSPOT (``RANGE first last`` / ``MULTIFRAME``).
    int range_first{0};
    int range_last{0};
    bool multiframe{false};
    /// ``≤ 0`` (default): keep per-frame spots for DPS; set ``MERGE_XY`` to cluster.
    double merge_xy_mm{0.0};
    double merge_phi_tol_deg{0.5};
    int min_frames{1};
    double phi_start_deg{0.0};
    double phi_end_deg{0.0};
    bool has_phi{false};
    float threshold_sigma{5.0f};
    int min_distance{8};
    /// Cap for indexing COLSPOT (large cells need ≫500 single-frame spots).
    int max_spots{1500};
    int border{20};
};

/// Keyword state accumulated by MOSFLM CONTROL (Phase D / E).
/// One AUTOINDEX lattice solution (a row of the penalty table).
struct AutoindexSolution {
    int rank{0};            ///< 1 = lowest penalty
    int character{0};       ///< Lattice character 1–44 (``AUTOINDEX CHARACTER n``)
    std::string lattice;    ///< Bravais lattice, e.g. ``"tP"``
    double penalty{0.0};    ///< Raw DPS goof (lower is better)
    UnitCell cell{};        ///< Conventional cell for that lattice
    std::optional<double> metric_deviation_pct;
    std::vector<std::string> space_groups;  ///< Possible space groups
};

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
    /// Default matches common fine-phi Eiger scans (≪ older 0.30°, which
    /// forced FRACTIONCALC≈0.5 on 0.15° frames when MOSAIC was omitted).
    double mosaic_deg{0.15};
    /// Beam divergence full-angles (degrees); MOSFLM ``DIVERGENCE`` / DIVH,DIVV.
    /// Default 0 unless ``DIVERGENCE`` is set (matches known-indexed insulin).
    double divh_deg{0.0};
    double divv_deg{0.0};
    double separation_x_mm{0.57};
    double separation_y_mm{0.57};
    /// High-resolution cutoff (Å). 0, or not ``dmin_explicit``, means PROCESS
    /// uses the farthest inset detector corner.
    double dmin_A{0.0};
    bool dmin_explicit{false};
    double dmax_A{200.0};
    /// Predictor φ sample step (degrees); 0.025 matches known_indexed_insulin.
    double phi_step_deg{0.025};
    /// Detector fiducials (neutral by default; insulin oracle via CONTROL
    /// ``YSCALE`` / ``TILT`` / ``TWIST`` / ``ROTMAT2``).
    double yscale{1.0};
    double tilt_mrad{0.0};
    double twist_mrad{0.0};
    bool use_rotmat2{false};
    /// If true, widen Ewald acceptance with mosaic (indexing-style). PROCESS
    /// default is false (REEK); set ``WIDEN_EWALD ON`` for insulin oracle.
    bool widen_ewald_with_mosaic{false};
    double mosaic_ewald_scale{0.55};
    /// Minimum FRACTIONCALC to keep on MTZ write (CONTROL ``MINFRAC``).
    /// Default 0.5 matches DIALS-style export; ≤ 0 disables.
    double min_fractioncalc{0.5};
    /// XDS ``FRACTION_OF_POLARIZATION`` for PROCESS LP (CONTROL
    /// ``POLARIZATION``). Default 0.98 (synchrotron); 0 = unpolarized.
    double polarization_fraction{0.98};
    /// XDS ``AIR=`` stand-in (mm⁻¹). Default small synchrotron air path.
    double air_absorption_mm_inv{0.0001};
    /// Scale partials by ``1/FRACTIONCALC`` (CONTROL ``PARTIALSCALE``).
    bool scale_by_partiality{true};
    /// Auto-estimate MOSAIC / DIVERGENCE from strong-spot widths.
    bool auto_sigma_msd{true};
    /// Rocking-curve prediction + partials summation in PROCESS (``ROCKING``).
    bool rocking_model{true};
    /// True once CONTROL ``MOSAIC`` was given (AUTOSIGMA then keeps it).
    bool mosaic_from_user{false};
    /// Auto-write GXPARM.XDS after PROCESS (default ON). Override path with
    /// ``XPARMOUT``; disable with ``XPARMOUT OFF``.
    bool auto_gxparm{true};
    /// Periodic post-refinement interval in frames (0 = only when
    /// ``POSTREFINE`` without PERIOD; CONTROL ``POSTREFINE PERIOD n``).
    int postrefine_period{25};
    /// Parameters positional refinement may change (CONTROL ``REFINE``).
    /// ``--xparm`` keeps the XDS-refined detector: ORIENTATION CELL only.
    bool refine_orientation{true};
    bool refine_cell{true};
    bool refine_distance{true};
    bool refine_beam{true};
    bool refine_tilt_twist{true};
    /// OpenMP worker count for MULTIFRAME FINDSPOTS / PROCESS (CONTROL
    /// ``NPROC``). Default from ``PROCESSINGLAB_NPROC`` or hardware.
    int nproc{mosflm::util::default_nproc()};

    UnitCell cell{78.19, 78.19, 78.19, 90.0, 90.0, 90.0};
    int symmetry_number{1};
    /// Apply MOSFLM-style CELLFIX before PROCESS/MTZ (default ON).
    bool cellfix_enabled{true};
    /// CELLFIX softness in (0,1]; 1 = hard metric projection.
    double cellfix_softness{1.0};
    /// MOSFLM Orientation Matrix [A] (lab X,Y,Z rows), before C++ frame map.
    Eigen::Matrix3d amat{Eigen::Matrix3d::Identity()};

    std::filesystem::path directory{"."};
    /// Where relative output paths go (spot list, XPARM.XDS, GXPARM.XDS);
    /// empty = ``directory``, as MOSFLM does. ``ProcessingLab run`` sets it
    /// so nothing is written next to the images.
    std::filesystem::path output_directory;
    std::string image_template;  ///< e.g. ``insulin_1_###.img``
    std::filesystem::path image_path;  ///< ``IMAGE`` keyword if set
    /// Output MTZ path (``HKLOUT`` keyword or CLI override).
    std::filesystem::path hklout{"process.mtz"};
    /// Optional GXPARM.XDS write path (``XPARMOUT``).
    std::filesystem::path xparm_out;
    bool has_xparm_out{false};

    bool process_pending{false};
    int process_first{1};
    int process_last{1};
    double process_start_phi_deg{0.0};
    double process_angle_deg{1.0};
    /// True when ``PROCESS … START`` / ``ANGLE`` overrode image φ metadata.
    bool has_process_phi{false};
    /// Added to every frame's φ start/end (CONTROL ``PHI_OFFSET``).
    double phi_offset_deg{0.0};
    /// If >= 0, overrides MEAS/profile peak half-width (CONTROL ``PEAKHALF``).
    int peak_half_px{-1};
    /// Peak half-width (px) measured from the strong spots at the start of
    /// PROCESS (MOSFLM/XDS-style); used when PEAKHALF is not set. 0 = none.
    int spot_peak_half_px{0};

    /// Standard-profile detector tiling (MOSFLM NMASKS / XLINES×YLINES).
    int nmasks_x{5};
    int nmasks_y{5};
    /// Opt-in Stage 2 post-refinement on PROCESS.
    bool enable_postrefinement{false};
    /// Shoebox depth in φ (odd). ``PROFILE 3D`` uses this as nφ; default 5
    /// when PROFILE 3D is set and NADD was left at 1. 2D PROCESS ignores it.
    int nadd{1};
    /// If true, PROCESS extracts (x,y,φ) shoeboxes and fits a 3D Rossmann
    /// profile (CONTROL ``PROFILE 3D``). Default remains 2D + ROCKING.
    bool profile_3d{false};

    bool findspots_pending{false};
    FindSpotsControlSettings findspots{};
    /// Number of spots written by the last FINDSPOTS GO.
    std::size_t findspots_n_spots{0};

    bool autoindex_pending{false};
    AutoIndexControlSettings autoindex{};
    /// Filled by AUTOINDEX runner for DiffractionLab / machine-readable dump.
    std::string autoindex_system;
    int autoindex_character{0};
    std::size_t autoindex_n_spots{0};
    /// Fraction of I/σ-filtered spots indexed by the chosen AMAT (0–1).
    double autoindex_indexed_frac{0.0};
    /// RMS fractional HKL residual for indexed spots.
    double autoindex_rmsd_hkl{0.0};
    /// Every lattice solution AUTOINDEX considered, best first.
    std::vector<AutoindexSolution> autoindex_solutions;
    /// Rank of the last solution before the largest penalty jump (LASTGOOD).
    int autoindex_lastgood_rank{0};
    /// INVERTX flag of the ``.spt`` AUTOINDEX indexed (header field 6).
    /// Selects the detector CCOM PROCESS must use so predictions land in the
    /// frame the AMAT was derived in (``ccom_from_spt_omega``). Default true
    /// = MOSFLM/ADSC convention, which is also what a user-supplied MATRIX
    /// assumes.
    bool spt_invertx{true};

    /// Multi-lattice: additional orientation matrices (first is ``amat``).
    std::vector<Eigen::Matrix3d> extra_lattices;

    bool strategy_pending{false};
    bool strategy_done{false};
    double strategy_osc_range_deg{90.0};
    double strategy_completeness_est{0.0};

    bool done_run{false};
    bool exit_requested{false};
};

/// ``name`` resolved for writing: as is if absolute, else under
/// ``output_directory`` (or ``directory`` when that is empty).
[[nodiscard]] inline std::filesystem::path output_path(
    const ControlState& state, const std::filesystem::path& name
) {
    if (name.is_absolute()) {
        return name;
    }
    return (state.output_directory.empty() ? state.directory : state.output_directory) / name;
}

/// Outcome of one ``GO`` / ``RUN``.
struct ControlGoResult {
    bool ran_process{false};
    bool ran_autoindex{false};
    bool ran_findspots{false};
    bool ran_strategy{false};
    diffraction::ProcessResult process{};
    crystal::CrystalMetadata crystal{};
    std::filesystem::path image_used{};
    /// Inclusive PROCESS image range actually integrated this GO.
    int process_first{0};
    int process_last{0};
    std::string message;
};

/// Hook that runs AUTOINDEX and writes CELL/MATRIX into ``state``.
/// ``image_override`` supplies header geometry CONTROL did not set.
using AutoIndexRunner = std::function<void(
    ControlState& state,
    const std::optional<std::filesystem::path>& image_override
)>;

/// Hook that runs FINDSPOTS (load image → peak find → write ``.spt``).
using FindSpotsRunner = std::function<void(
    ControlState& state,
    const std::optional<std::filesystem::path>& image_override
)>;

/// Hooks around each ``GO`` in :meth:`ControlSession::run_stream`, for run
/// records: ``before`` sees the state with the step still pending, ``after``
/// the result. An exception from ``go`` propagates after ``before`` ran.
struct GoObserver {
    std::function<void(const ControlState&)> before;
    std::function<void(const ControlGoResult&, const ControlState&)> after;
};

/// Pure-C++ CONTROL driver (Phase D subset + Phase E extensions).
///
/// Supported: geometry keywords, ``PROCESS``, ``GO``/``RUN``, ``EXIT``,
/// ``FINDSPOTS`` / ``AUTOINDEX`` (when runners are registered), ``STRATEGY``,
/// ``NMASKS``, ``POSTREFINE``, ``NADD``, ``MULTILATTICE``. Unimplemented major
/// verbs (``REFIX``, …) throw. Other unknown keywords warn once then continue.
class ControlSession {
public:
    /// Register DPS AUTOINDEX implementation (from ``processinglab_indexing``).
    static void register_autoindex_runner(AutoIndexRunner runner);

    /// Register DiffractionLab-style FINDSPOTS (from ``processinglab_indexing``).
    static void register_findspots_runner(FindSpotsRunner runner);

    /// Parse one logical line (comments stripped). Throws on malformed
    /// known keywords or unimplemented major verbs.
    void parse_line(std::string_view line);

    /// Parse a whole script file.
    void parse_file(const std::filesystem::path& path);

    /// Execute ``GO``/``RUN``: AUTOINDEX, STRATEGY, and/or PROCESS.
    [[nodiscard]] ControlGoResult go(
        const std::optional<std::filesystem::path>& image_override = std::nullopt
    );

    /// Build a PROCESS request for one 1-based image number (φ / batch).
    [[nodiscard]] diffraction::ProcessRequest make_process_request(
        const image::DiffractionImage& image,
        int frame
    ) const;

    /// Resolve the path for ``frame`` (TEMPLATE expands; HDF5 masters reuse path).
    [[nodiscard]] std::filesystem::path resolve_image_path(
        const std::optional<std::filesystem::path>& image_override =
            std::nullopt,
        int frame = -1
    ) const;

    [[nodiscard]] const ControlState& state() const { return state_; }
    /// Set :member:`ControlState::output_directory`.
    void set_output_directory(const std::filesystem::path& dir) { state_.output_directory = dir; }
    /// Observe each ``GO`` run by :meth:`run_stream`.
    void set_go_observer(GoObserver observer) { go_observer_ = std::move(observer); }
    [[nodiscard]] bool exit_requested() const { return state_.exit_requested; }
    [[nodiscard]] bool done_run() const { return state_.done_run; }
    [[nodiscard]] bool strategy_done() const { return state_.strategy_done; }

    /// True when EOF should emit Fortran-style no-RUN warning.
    [[nodiscard]] bool needs_no_run_warning() const {
        return !state_.done_run && !state_.strategy_done;
    }

    /// Read CONTROL lines from ``in`` until EOF or EXIT, running each GO/RUN.
    ///
    /// \param in             Command source (stdin, script file, …).
    /// \param image_override Image used instead of IMAGE/TEMPLATE if set.
    /// \param prompt         If non-empty, written to stdout and flushed
    ///     before each line is read (interactive terminals only).
    [[nodiscard]] std::vector<ControlGoResult> run_stream(
        std::istream& in,
        const std::optional<std::filesystem::path>& image_override =
            std::nullopt,
        std::string_view prompt = {}
    );

    [[nodiscard]] std::vector<ControlGoResult> run_script(
        const std::filesystem::path& script,
        const std::optional<std::filesystem::path>& image_override =
            std::nullopt
    );

private:
    ControlState state_;
    GoObserver go_observer_;
    int matrix_needed_{0};
    std::vector<double> matrix_buf_;
    std::set<std::string> warned_unknown_;

    void apply_keyword(const std::string& key, const std::vector<std::string>& args);
    void ingest_matrix_numbers(const std::vector<double>& vals);
    void finish_matrix_if_ready();
    void parse_autoindex_args(const std::vector<std::string>& args);
    void parse_findspots_args(const std::vector<std::string>& args);
};

}  // namespace mosflm::control
