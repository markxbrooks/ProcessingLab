#pragma once

#include "mosflm/io/ExperimentModel.hpp"

#include <filesystem>
#include <string>

namespace mosflm::io {

/// Files written at the start of a run so it can be repeated in another run
/// mode (or by XDS / xia2 / a CONTROL script), each in its own subdirectory.
struct RepeatInputs {
    std::filesystem::path xds_inp;        ///< ``xds/XDS.INP``
    std::filesystem::path xinfo;          ///< ``xia2/automatic.xinfo``
    std::filesystem::path xia2_script;    ///< ``xia2/run_xia2.sh``
    std::filesystem::path mosflm_script;  ///< ``mosflm/mosflm.in``
};

/// ``XDS.INP`` for the resolved experiment: images, frame and spot ranges,
/// detector, geometry and any known cell / space group. ProcessingLab reads
/// it back with ``--run-mode xds``; for XDS itself the detector keywords
/// (``DETECTOR``, ``OVERLOAD``, ``ROTATION_AXIS``) are derived from the header
/// and should be checked.
///
/// :param m: The resolved experiment.
/// :return: The file text.
[[nodiscard]] std::string xds_inp_text(const ExperimentModel& m);

/// Write ``xds_inp_text(model)`` to ``path`` (atomically).
///
/// :param path: Destination, usually ``<output>/xds/XDS.INP``.
/// :param model: The experiment after every source has been merged.
/// :raises std::runtime_error: If the file cannot be written.
void write_xds_inp(const std::filesystem::path& path, const ExperimentModel& model);

/// xia2 ``.xinfo`` for the resolved experiment: one crystal and sweep, with
/// ``USER_CELL`` / ``USER_SPACEGROUP`` when known and ``BEAM`` / ``DISTANCE``
/// only where the run overrode the image header.
///
/// :param m: The resolved experiment.
/// :return: The file text.
[[nodiscard]] std::string xinfo_text(const ExperimentModel& m);

/// CONTROL script that repeats the run with ``--run-mode mosflm``: the
/// run's own CONTROL, self-contained (``IMAGE`` for an HDF5 master) and
/// without ``HKLOUT``, so outputs go to the new run's ``--output-dir``.
///
/// :param m: The resolved experiment.
/// :param control: The CONTROL text the run executes.
/// :return: The file text.
[[nodiscard]] std::string mosflm_script_text(const ExperimentModel& m, const std::string& control);

/// Write ``xds/``, ``xia2/`` and ``mosflm/`` under ``dir``.
///
/// :raises std::runtime_error: If a file cannot be written.
RepeatInputs write_repeat_inputs(
    const std::filesystem::path& dir, const ExperimentModel& m, const std::string& control
);

/// After the run: copy ``XPARM.XDS`` and ``GXPARM.XDS`` (those that exist)
/// from ``dir`` into ``dir/xds``.
void copy_repeat_geometry(const std::filesystem::path& dir);

}  // namespace mosflm::io
