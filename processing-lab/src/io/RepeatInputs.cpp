#include "mosflm/io/RepeatInputs.hpp"

#include "mosflm/crystal/SpaceGroupInfo.hpp"
#include "mosflm/diffraction/Defpix.hpp"
#include "mosflm/io/Json.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <system_error>

namespace mosflm::io {
namespace {

std::string num(double v, int prec) {
    std::ostringstream s;
    s << std::fixed << std::setprecision(prec) << v;
    return s.str();
}

/// ``v`` to full double precision, so a repeat run reads back the same
/// number: the beam centre feeds IDXREF, whose search over a flat optimum
/// (one frame, few spots) ends elsewhere for a 1e-6 mm difference.
std::string exact(double v) {
    std::ostringstream s;
    s << std::setprecision(std::numeric_limits<double>::max_digits10) << v;
    return s.str();
}

std::string upper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return s;
}

/// XDS ``DETECTOR=`` name and a default overload, from the header id.
std::pair<std::string, double> xds_detector(const ExperimentModel& m) {
    const std::string id = upper(m.detector_id);
    if (id.find("EIGER") != std::string::npos) return {"EIGER", 1048576.0};
    if (id.find("PILATUS") != std::string::npos) return {"PILATUS", 1048500.0};
    if (id.find("ADSC") != std::string::npos) return {"ADSC", 65535.0};
    if (id.find("RAYONIX") != std::string::npos || id.find("MAR") != std::string::npos) {
        return {"CCDCHESS", 65535.0};
    }
    return {"", 0.0};
}

/// ``ROTATION_AXIS`` from the header's axis omega (270: along +fast).
std::string xds_rotation_axis(const ExperimentModel& m) {
    if (!m.hdf5 && m.image_template.find(".img") != std::string::npos) {
        return "1 0 0";  // ADSC SMV: XDS reads the frame transposed
    }
    const double w = std::fmod(std::fmod(m.rotation_axis_omega_deg, 360.0) + 360.0, 360.0);
    if (std::fabs(w - 90.0) < 45.0) return "-1 0 0";
    if (std::fabs(w - 180.0) < 45.0) return "0 -1 0";
    if (std::fabs(w - 270.0) < 45.0) return "1 0 0";
    return "0 1 0";
}

/// ``x_master.h5`` → ``x_??????.h5`` (XDS names Eiger data files).
std::filesystem::path xds_template(const ExperimentModel& m) {
    if (m.hdf5) {
        const std::string name = m.image.value.filename().string();
        const std::string suffix = "master.h5";
        if (name.size() > suffix.size() && name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0) {
            return m.image.value.parent_path() / (name.substr(0, name.size() - suffix.size()) + "??????.h5");
        }
        return m.image.value;
    }
    std::string t = m.image_template;
    std::replace(t.begin(), t.end(), '#', '?');
    return m.directory / t;
}

/// ``P 41 21 2`` → ``P41212`` (xia2 symbols have no spaces).
std::string compact_symbol(int number) {
    const auto* sg = crystal::SpaceGroupInfo::find_by_number(number);
    if (sg == nullptr) {
        return {};
    }
    std::string s = sg->hm;
    s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
    return s;
}

void write_text(const std::filesystem::path& path, const std::string& text) {
    std::filesystem::create_directories(path.parent_path());
    write_file_atomic(path, text);
}

}  // namespace

std::string xds_inp_text(const ExperimentModel& m) {
    const auto [f0, f1] = m.frame_range.value;
    const auto [detector, default_overload] = xds_detector(m);
    std::ostringstream o;
    o << "! XDS.INP written by ProcessingLab at the start of the run.\n"
      << "! Repeat in ProcessingLab:  ProcessingLab run --run-mode xds --xds-inp XDS.INP\n"
      << "! For XDS, check DETECTOR, OVERLOAD, ROTATION_AXIS and the trusted region.\n"
      << "JOB= XYCORR INIT COLSPOT IDXREF DEFPIX INTEGRATE CORRECT\n";
    if (!detector.empty()) {
        const double overload = m.saturation_counts > 0.0 ? m.saturation_counts : default_overload;
        o << "DETECTOR=" << detector << " MINIMUM_VALID_PIXEL_VALUE=0 OVERLOAD="
          << static_cast<long long>(overload) << '\n';
    } else {
        o << "! DETECTOR= (not recognised from header '" << m.detector_id << "')\n";
    }
    o << "NAME_TEMPLATE_OF_DATA_FRAMES=" << xds_template(m).string() << '\n';
    if (m.hdf5) {
        o << "! LIB= /path/to/durin-plugin.so   (XDS reads HDF5 through a plugin)\n";
    }
    o << "DATA_RANGE=" << f0 << ' ' << f1 << '\n'
      << "SPOT_RANGE=" << f0 << ' ' << std::min(f0 + 49, f1) << '\n'
      << "BACKGROUND_RANGE=" << f0 << ' ' << std::min(f0 + 4, f1) << '\n';
    if (m.oscillation_deg.set()) {
        o << "OSCILLATION_RANGE=" << num(m.oscillation_deg.value, 4) << '\n';
    }
    o << "STARTING_ANGLE=" << num(m.phi_start_deg, 4) << " STARTING_FRAME=" << (m.hdf5 ? 1 : m.first_frame)
      << '\n';
    if (m.wavelength_A.set()) {
        o << "X-RAY_WAVELENGTH=" << num(m.wavelength_A.value, 6) << '\n';
    }
    if (m.distance_mm.set()) {
        o << "DETECTOR_DISTANCE=" << num(m.distance_mm.value, 3) << '\n';
    }
    if (m.pixel_mm.set()) {
        o << "NX=" << m.detector_nx << " NY=" << m.detector_ny << " QX=" << num(m.pixel_mm.value.first, 6)
          << " QY=" << num(m.pixel_mm.value.second, 6) << '\n';
        if (m.beam_mm.set()) {
            o << "ORGX=" << exact(m.beam_mm.value.first / m.pixel_mm.value.first)
              << " ORGY=" << exact(m.beam_mm.value.second / m.pixel_mm.value.second) << '\n';
        }
    }
    o << "ROTATION_AXIS=" << xds_rotation_axis(m) << '\n'
      << "DIRECTION_OF_DETECTOR_X-AXIS=1 0 0\n"
      << "DIRECTION_OF_DETECTOR_Y-AXIS=0 1 0\n"
      << "INCIDENT_BEAM_DIRECTION=0 0 1\n"
      << "FRACTION_OF_POLARIZATION=0.99\n"
      << "POLARIZATION_PLANE_NORMAL=0 1 0\n"
      << "FRIEDEL'S_LAW=TRUE\n";
    if (m.symmetry.set()) {
        o << "SPACE_GROUP_NUMBER=" << m.symmetry.value << '\n';
    } else {
        o << "! SPACE_GROUP_NUMBER= (from autoindexing)\n";
    }
    if (m.cell.set()) {
        o << "UNIT_CELL_CONSTANTS=";
        for (const double v : m.cell.value) {
            o << ' ' << num(v, 3);
        }
        o << '\n';
    }
    // The run's untrusted regions and excluded resolution (CONTROL → XDS).
    for (const auto& line : m.imported.trusted_region_control) {
        std::istringstream s(line);
        std::string key;
        std::string what;
        s >> key >> what;
        std::string rest;
        std::getline(s, rest);
        if (key == "UNTRUSTED") {
            o << "UNTRUSTED_" << what << '=' << rest << '\n';
        } else if (key == "EXCLUDE" && what == "ICE") {
            for (const auto& r : diffraction::ice_ring_ranges()) {
                o << "EXCLUDE_RESOLUTION_RANGE= " << num(r.d1_A, 3) << ' ' << num(r.d2_A, 3) << '\n';
            }
        } else if (key == "EXCLUDE") {
            o << "EXCLUDE_RESOLUTION_RANGE=" << rest << '\n';
        } else if (key == "DEFPIX" && what == "VALUE_RANGE") {
            o << "VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS=" << rest << '\n';
        }
    }
    return o.str();
}

void write_xds_inp(const std::filesystem::path& path, const ExperimentModel& model) {
    write_text(path, xds_inp_text(model));
}

std::string xinfo_text(const ExperimentModel& m) {
    const auto [f0, f1] = m.frame_range.value;
    const auto overridden = [](MetadataSource s) {
        return s != MetadataSource::None && s != MetadataSource::ImageHeader;
    };
    std::ostringstream o;
    o << "! xia2 .xinfo written by ProcessingLab at the start of the run.\n"
      << "! Repeat in ProcessingLab:  ProcessingLab run --run-mode xia2 --xinfo automatic.xinfo\n"
      << "! BEAM / DISTANCE appear only where the run overrode the image header\n"
      << "! (BEAM in ProcessingLab's frame, ORGX×QX ORGY×QY: check before xia2 uses it).\n"
      << "BEGIN PROJECT PROCESSINGLAB\n"
      << "BEGIN CRYSTAL XTAL\n";
    if (m.symmetry.set()) {
        if (const std::string sym = compact_symbol(m.symmetry.value); !sym.empty()) {
            o << "USER_SPACEGROUP " << sym << '\n';
        }
    }
    if (m.cell.set()) {
        o << "USER_CELL";
        for (const double v : m.cell.value) {
            o << ' ' << num(v, 3);
        }
        o << '\n';
    }
    o << "BEGIN WAVELENGTH NATIVE\n";
    if (m.wavelength_A.set()) {
        o << "WAVELENGTH " << num(m.wavelength_A.value, 6) << '\n';
    }
    o << "END WAVELENGTH NATIVE\n"
      << "BEGIN SWEEP SWEEP1\n"
      << "WAVELENGTH NATIVE\n"
      << "DIRECTORY " << m.image.value.parent_path().string() << '\n'
      << "IMAGE " << m.image.value.filename().string() << '\n'
      << "START_END " << f0 << ' ' << f1 << '\n';
    if (m.beam_mm.set() && overridden(m.beam_mm.source)) {
        o << "BEAM " << exact(m.beam_mm.value.first) << ' ' << exact(m.beam_mm.value.second) << '\n';
    }
    if (m.distance_mm.set() && overridden(m.distance_mm.source)) {
        o << "DISTANCE " << num(m.distance_mm.value, 3) << '\n';
    }
    o << "END SWEEP SWEEP1\n"
      << "END CRYSTAL XTAL\n"
      << "END PROJECT PROCESSINGLAB\n";
    return o.str();
}

std::string mosflm_script_text(const ExperimentModel& m, const std::string& control) {
    std::ostringstream o;
    o << "! CONTROL script written by ProcessingLab at the start of the run.\n"
      << "! Repeat:  ProcessingLab --run-mode mosflm --output-dir NEW_DIR mosflm.in\n"
      << "! MOSFLM keywords plus ProcessingLab's FINDSPOTS; outputs go to --output-dir.\n";
    if (m.hdf5) {
        o << "IMAGE " << m.image.value.string() << '\n';
    }
    std::istringstream in(control);
    std::string line;
    while (std::getline(in, line)) {
        std::string key = line.substr(0, line.find(' '));
        if (upper(key) == "HKLOUT") {
            continue;  // the new run's --output-dir decides where the MTZ goes
        }
        o << line << '\n';
    }
    return o.str();
}

RepeatInputs write_repeat_inputs(
    const std::filesystem::path& dir, const ExperimentModel& m, const std::string& control
) {
    RepeatInputs r{dir / "xds" / "XDS.INP", dir / "xia2" / "automatic.xinfo", dir / "xia2" / "run_xia2.sh",
                   dir / "mosflm" / "mosflm.in"};
    write_xds_inp(r.xds_inp, m);
    write_text(r.xinfo, xinfo_text(m));
    write_text(r.xia2_script,
               "#!/bin/sh\n"
               "# Repeat this run with xia2 (written by ProcessingLab).\n"
               "# In ProcessingLab:  ProcessingLab run --run-mode xia2 --xinfo automatic.xinfo\n"
               "cd \"$(dirname \"$0\")\" && exec xia2 xinfo=automatic.xinfo \"$@\"\n");
    std::error_code ec;
    std::filesystem::permissions(r.xia2_script,
                                 std::filesystem::perms::owner_exec | std::filesystem::perms::group_exec |
                                     std::filesystem::perms::others_exec,
                                 std::filesystem::perm_options::add, ec);
    write_text(r.mosflm_script, mosflm_script_text(m, control));
    return r;
}

void copy_repeat_geometry(const std::filesystem::path& dir) {
    for (const char* name : {"XPARM.XDS", "GXPARM.XDS"}) {
        std::error_code ec;
        if (std::filesystem::exists(dir / name, ec)) {
            std::filesystem::create_directories(dir / "xds", ec);
            std::filesystem::copy_file(dir / name, dir / "xds" / name,
                                       std::filesystem::copy_options::overwrite_existing, ec);
        }
    }
}

}  // namespace mosflm::io
