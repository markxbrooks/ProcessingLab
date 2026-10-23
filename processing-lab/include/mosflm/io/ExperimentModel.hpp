#pragma once

#include "mosflm/io/ImportedSweepGeometry.hpp"
#include "mosflm/io/Json.hpp"

#include <array>
#include <ctime>
#include <filesystem>
#include <optional>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

namespace mosflm::io {

/// Where a resolved experiment value came from, lowest priority first.
enum class MetadataSource {
    None,
    ImageHeader,
    XdsInp,
    Xinfo,
    Xparm,    ///< ``--xparm`` XPARM.XDS / GXPARM.XDS
    PreviousRun,  ///< ``run --from``: the earlier run's refined GXPARM.XDS
    Control,  ///< CONTROL script keyword (script mode)
    CommandLine,
};

/// Label for logs (``"HDF5 header"``, ``"XDS.INP"``, …).
[[nodiscard]] const char* source_label(MetadataSource source, bool hdf5);
/// Identifier for records (``"image_header"``, ``"xds_inp"``, …).
[[nodiscard]] const char* source_id(MetadataSource source);

/// A value with the source that set it, and what it replaced.
template <class T>
struct Sourced {
    T value{};
    MetadataSource source{MetadataSource::None};
    /// Earlier values, oldest first, each with its source: a later source
    /// replaces the value but the record keeps what it overrode.
    std::vector<std::pair<T, MetadataSource>> overrides;

    [[nodiscard]] bool set() const { return source != MetadataSource::None; }
    /// Take ``v`` from ``src``, keeping the previous value in ``overrides``.
    void assign(T v, MetadataSource src) {
        if (set()) {
            overrides.emplace_back(value, source);
        }
        value = std::move(v);
        source = src;
    }
};

/// Inputs for one sweep of ``ProcessingLab run``.
struct ExperimentInputs {
    /// ``--image``: an HDF5/NeXus master, or one frame of a CBF/SMV sweep.
    std::optional<std::filesystem::path> image;
    std::optional<std::filesystem::path> xds_inp;
    std::optional<std::filesystem::path> xinfo;
    /// ``--sweep``: which ``.xinfo`` sweep (required if it has several).
    std::optional<std::string> sweep;
    /// ``--range FIRST:LAST`` (inclusive).
    std::optional<std::pair<int, int>> range;
    /// ``--xparm`` (or the GXPARM.XDS of ``--from``): cell, symmetry,
    /// geometry and orientation, above ``.xinfo`` and below the command line.
    std::optional<std::filesystem::path> xparm;
    /// Provenance of ``xparm``: ``Xparm``, or ``PreviousRun`` for ``--from``.
    MetadataSource xparm_source{MetadataSource::Xparm};
    /// Detector tilt and twist (mrad) from ``--from``'s ``process.json``
    /// (XPARM files do not carry them); provenance ``xparm_source``.
    std::optional<std::pair<double, double>> tilt_twist_mrad;
};

/// The experiment once every source is merged, with per-field provenance.
///
/// Precedence per field, lowest to highest: image header, ``XDS.INP``,
/// ``.xinfo``, XPARM (``--xparm`` / ``--from``), command line. A source only overrides the fields it sets, so
/// an ``.xinfo`` without detector geometry keeps the header's.
struct ExperimentModel {
    Sourced<std::filesystem::path> image;  ///< HDF5 master or first frame.
    bool hdf5{false};                      ///< ``image`` holds every frame.
    std::filesystem::path directory;       ///< CBF/SMV: frame directory.
    std::string image_template;            ///< CBF/SMV: ``name_####.cbf``.
    Sourced<double> wavelength_A;
    Sourced<double> distance_mm;
    Sourced<std::pair<double, double>> beam_mm;
    Sourced<std::pair<double, double>> pixel_mm;
    Sourced<double> oscillation_deg;
    Sourced<std::array<double, 6>> cell;
    Sourced<int> symmetry;
    Sourced<std::pair<int, int>> frame_range;
    /// Orientation matrix (row-major) from an XPARM; when set, ``run`` skips
    /// autoindexing and refines only the orientation.
    Sourced<std::array<double, 9>> orientation;
    /// Detector tilt / twist (mrad), when an earlier run refined them.
    Sourced<std::pair<double, double>> tilt_twist_mrad;
    int frames_available{0};  ///< Frames found on disk / in the master.
    /// Scientifically meaningful disagreements between sources (the higher
    /// priority value is used; nothing is reconciled silently).
    std::vector<std::string> warnings;
    /// Values from imports/CLI, for CONTROL (header values are left to the
    /// session's own header defaults).
    ImportedSweepGeometry imported;
    /// CONTROL lines that set ``orientation`` (``INVERTX``, ``REFINE``,
    /// ``MATRIX``), after the imported geometry.
    std::vector<std::string> orientation_control;

    /// CONTROL lines: DIRECTORY/TEMPLATE for CBF/SMV, then imported geometry.
    [[nodiscard]] std::vector<std::string> control_lines() const;
    /// "Experiment metadata" block for the log.
    void write_summary(std::ostream& os) const;
    /// ``experiment`` section of ``input.json``: each field as
    /// ``{value, unit, source, overrides: [{value, source}, …]}``.
    [[nodiscard]] Json to_json() const;
};

/// Script mode: geometry set in a parsed (not run) CONTROL ``state`` that
/// differs from ``model`` — from ``--xparm`` or the script's keywords —
/// overrides it with provenance ``source``; differences are warned like
/// imports. Values the imports already supplied are unchanged.
void apply_control_state(
    ExperimentModel& model,
    const control::ControlState& state,
    MetadataSource source = MetadataSource::Control
);

/// Merge image header, ``XDS.INP``, ``.xinfo`` and CLI values for one sweep.
///
/// :param inputs: Command-line inputs (``image`` may come from an import).
/// :return: The resolved model.
/// :raises std::runtime_error: If no image source is given, the ``.xinfo``
///     sweep is ambiguous or missing, or the image cannot be read.
[[nodiscard]] ExperimentModel resolve_experiment(const ExperimentInputs& inputs);

/// Write ``path`` as an XDS ``XDS.INP`` of the resolved options, at the
/// start of a run. Only fields that are set are written, so a later
/// ``--xds-inp`` of this file reproduces them. CBF/SMV ``#`` placeholders
/// become ``?``. An HDF5 master is named as ``NAME_TEMPLATE_OF_DATA_FRAMES``.
/// Orientation is omitted (that belongs in ``XPARM.XDS``).
///
/// :param path: Destination, usually ``<output>/XDS.INP``.
/// :param model: The experiment after every source has been merged.
/// :raises std::runtime_error: If the file cannot be written.
void write_xds_inp(const std::filesystem::path& path, const ExperimentModel& model);

/// Create a fresh output directory and return it.
///
/// With ``requested`` it must not exist yet (no silent reuse). Otherwise
/// ``parent/processing_YYYYMMDD_HHMMSS`` is created, with ``_2``, ``_3``, …
/// appended if another run took that name first (creation is atomic, so two
/// simultaneous jobs never share a directory).
///
/// :raises std::runtime_error: If ``requested`` exists or nothing can be made.
[[nodiscard]] std::filesystem::path create_output_directory(
    const std::optional<std::filesystem::path>& requested,
    const std::filesystem::path& parent,
    std::time_t now
);

}  // namespace mosflm::io
