#pragma once

#include "mosflm/io/ImportedSweepGeometry.hpp"

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
    CommandLine,
};

/// Label for logs (``"HDF5 header"``, ``"XDS.INP"``, …).
[[nodiscard]] const char* source_label(MetadataSource source, bool hdf5);

/// A value with the source that set it.
template <class T>
struct Sourced {
    T value{};
    MetadataSource source{MetadataSource::None};

    [[nodiscard]] bool set() const { return source != MetadataSource::None; }
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
};

/// The experiment once every source is merged, with per-field provenance.
///
/// Precedence per field, lowest to highest: image header, ``XDS.INP``,
/// ``.xinfo``, command line. A source only overrides the fields it sets, so
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
    int frames_available{0};  ///< Frames found on disk / in the master.
    /// Scientifically meaningful disagreements between sources (the higher
    /// priority value is used; nothing is reconciled silently).
    std::vector<std::string> warnings;
    /// Values from imports/CLI, for CONTROL (header values are left to the
    /// session's own header defaults).
    ImportedSweepGeometry imported;

    /// CONTROL lines: DIRECTORY/TEMPLATE for CBF/SMV, then imported geometry.
    [[nodiscard]] std::vector<std::string> control_lines() const;
    /// "Experiment metadata" block for the log.
    void write_summary(std::ostream& os) const;
};

/// Merge image header, ``XDS.INP``, ``.xinfo`` and CLI values for one sweep.
///
/// :param inputs: Command-line inputs (``image`` may come from an import).
/// :return: The resolved model.
/// :raises std::runtime_error: If no image source is given, the ``.xinfo``
///     sweep is ambiguous or missing, or the image cannot be read.
[[nodiscard]] ExperimentModel resolve_experiment(const ExperimentInputs& inputs);

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
