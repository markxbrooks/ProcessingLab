#include "mosflm/io/ExperimentModel.hpp"

#include "mosflm/image/ImageIO.hpp"
#include "mosflm/io/XdsInpImporter.hpp"
#include "mosflm/io/XinfoImporter.hpp"
#include "mosflm/io/XparmIO.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <system_error>

namespace mosflm::io {
namespace {

bool is_container(const std::filesystem::path& p) {
    std::string ext = p.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return ext == ".h5" || ext == ".hdf5" || ext == ".nxs";
}

/// Frame number in a CBF/SMV file name (last digit run of the stem), or 1.
int frame_number(const std::filesystem::path& p) {
    const std::string stem = p.stem().string();
    const auto end = stem.find_last_of("0123456789");
    if (end == std::string::npos) {
        return 1;
    }
    auto begin = end;
    while (begin > 0 && std::isdigit(static_cast<unsigned char>(stem[begin - 1]))) {
        --begin;
    }
    return std::stoi(stem.substr(begin, end - begin + 1));
}

/// ``name_####.cbf`` with ``frame`` substituted.
std::filesystem::path frame_path(
    const std::filesystem::path& dir, const std::string& tmpl, int frame
) {
    const auto hash = tmpl.find('#');
    if (hash == std::string::npos) {
        return dir / tmpl;
    }
    std::size_t n = 0;
    while (hash + n < tmpl.size() && tmpl[hash + n] == '#') {
        ++n;
    }
    std::string num = std::to_string(frame);
    if (num.size() < n) {
        num.insert(0, n - num.size(), '0');
    }
    return dir / (tmpl.substr(0, hash) + num + tmpl.substr(hash + n));
}

std::string fmt(double v, int prec) {
    std::ostringstream o;
    o << std::fixed << std::setprecision(prec) << v;
    return o.str();
}

/// Overlay ``value`` from ``source``; warn if it changes an earlier source's
/// value by more than ``tol`` (relative if ``rel``).
template <class T, class Diff>
void overlay(
    Sourced<T>& field, const T& value, MetadataSource source, bool hdf5,
    const char* what, Diff&& differs, std::vector<std::string>& warnings,
    const std::string& shown_new, const std::string& shown_old
) {
    if (field.set() && differs(field.value, value)) {
        warnings.push_back(
            std::string(what) + ": " + source_label(source, hdf5) + " " + shown_new +
            " overrides " + source_label(field.source, hdf5) + " " + shown_old
        );
    }
    field.assign(value, source);
}

std::string cell_str(const std::array<double, 6>& c) {
    std::string out;
    for (double v : c) {
        out += (out.empty() ? "" : " ") + fmt(v, 2);
    }
    return out;
}

bool rel_differs(double a, double b, double tol) {
    return std::fabs(a - b) > tol * std::max(std::fabs(a), std::fabs(b));
}

/// Fold one imported source into the model, field by field.
void apply_import(
    ExperimentModel& m, const ImportedSweepGeometry& g, MetadataSource src
) {
    auto& w = m.warnings;
    const bool h = m.hdf5;
    if (g.has_wavelength) {
        overlay(m.wavelength_A, g.wavelength_A, src, h, "Wavelength",
                [](double a, double b) { return rel_differs(a, b, 1e-3); }, w,
                fmt(g.wavelength_A, 5) + " Å", fmt(m.wavelength_A.value, 5) + " Å");
    }
    if (g.has_distance) {
        overlay(m.distance_mm, g.distance_mm, src, h, "Distance",
                [](double a, double b) { return rel_differs(a, b, 5e-3); }, w,
                fmt(g.distance_mm, 2) + " mm", fmt(m.distance_mm.value, 2) + " mm");
    }
    if (g.has_beam) {
        const std::pair<double, double> v{g.beam_x_mm, g.beam_y_mm};
        overlay(m.beam_mm, v, src, h, "Beam",
                [](const auto& a, const auto& b) {
                    return std::hypot(a.first - b.first, a.second - b.second) > 0.5;
                }, w,
                fmt(v.first, 2) + " " + fmt(v.second, 2) + " mm",
                fmt(m.beam_mm.value.first, 2) + " " + fmt(m.beam_mm.value.second, 2) + " mm");
    }
    if (g.has_pixel) {
        const std::pair<double, double> v{g.pixel_x_mm, g.pixel_y_mm};
        overlay(m.pixel_mm, v, src, h, "Pixel",
                [](const auto& a, const auto& b) {
                    return rel_differs(a.first, b.first, 1e-2) || rel_differs(a.second, b.second, 1e-2);
                }, w,
                fmt(v.first, 4) + " " + fmt(v.second, 4) + " mm",
                fmt(m.pixel_mm.value.first, 4) + " " + fmt(m.pixel_mm.value.second, 4) + " mm");
    }
    if (g.has_oscillation && g.oscillation_deg > 0.0) {
        overlay(m.oscillation_deg, g.oscillation_deg, src, h, "Oscillation",
                [](double a, double b) { return rel_differs(a, b, 1e-2); }, w,
                fmt(g.oscillation_deg, 4) + "°", fmt(m.oscillation_deg.value, 4) + "°");
    }
    if (g.has_cell) {
        const std::array<double, 6> v{g.cell_a, g.cell_b, g.cell_c,
                                      g.cell_alpha, g.cell_beta, g.cell_gamma};
        overlay(m.cell, v, src, h, "Cell",
                [](const auto& a, const auto& b) {
                    for (int k = 0; k < 6; ++k) {
                        if (rel_differs(a[k], b[k], 1e-2)) {
                            return true;
                        }
                    }
                    return false;
                }, w, cell_str(v), cell_str(m.cell.value));
    }
    if (g.has_symmetry) {
        overlay(m.symmetry, g.symmetry_number, src, h, "Space group",
                [](int a, int b) { return a != b; }, w,
                std::to_string(g.symmetry_number), std::to_string(m.symmetry.value));
    }
    if (g.has_process_range) {
        overlay(m.frame_range, std::pair<int, int>{g.process_first, g.process_last},
                src, h, "Frame range", [](const auto&, const auto&) { return false; },
                w, "", "");
    }
    m.imported.merge_from(g);
}

/// Image source from an import (directory + image or template).
std::optional<std::filesystem::path> image_from_import(const ImportedSweepGeometry& g) {
    const std::filesystem::path dir = g.has_directory ? g.directory : std::filesystem::path(".");
    if (g.has_image) {
        return g.image_path.is_absolute() ? g.image_path : dir / g.image_path;
    }
    if (g.has_template) {
        const std::string t = xds_question_marks_to_hash(g.image_template);
        return frame_path(dir, t, g.has_process_range ? g.process_first : 1);
    }
    return std::nullopt;
}

}  // namespace

const char* source_label(MetadataSource source, bool hdf5) {
    switch (source) {
        case MetadataSource::ImageHeader:
            return hdf5 ? "HDF5 header" : "image header";
        case MetadataSource::XdsInp:
            return "XDS.INP";
        case MetadataSource::Xinfo:
            return ".xinfo";
        case MetadataSource::Xparm:
            return "XPARM";
        case MetadataSource::PreviousRun:
            return "previous run";
        case MetadataSource::Control:
            return "CONTROL";
        case MetadataSource::CommandLine:
            return "CLI";
        case MetadataSource::None:
            break;
    }
    return "unset";
}

const char* source_id(MetadataSource source) {
    switch (source) {
        case MetadataSource::ImageHeader: return "image_header";
        case MetadataSource::XdsInp: return "xds_inp";
        case MetadataSource::Xinfo: return "xinfo";
        case MetadataSource::Xparm: return "xparm";
        case MetadataSource::PreviousRun: return "previous_run";
        case MetadataSource::Control: return "control";
        case MetadataSource::CommandLine: return "command_line";
        case MetadataSource::None: break;
    }
    return "unset";
}

ExperimentModel resolve_experiment(const ExperimentInputs& in) {
    ExperimentModel m;
    std::optional<ImportedSweepGeometry> xds;
    std::optional<ImportedSweepGeometry> xinfo;
    if (in.xds_inp) {
        xds = XdsInpImporter::parse_file(*in.xds_inp);
    }
    if (in.xinfo) {
        const auto names = XinfoImporter::sweep_names_in_file(*in.xinfo);
        if (in.sweep) {
            xinfo = XinfoImporter::parse_file(*in.xinfo, *in.sweep);
        } else if (names.size() > 1) {
            std::string list;
            for (const auto& n : names) {
                list += (list.empty() ? "" : ", ") + n;
            }
            throw std::runtime_error(
                in.xinfo->string() + " has " + std::to_string(names.size()) +
                " sweeps (" + list + "); choose one with --sweep"
            );
        } else {
            xinfo = XinfoImporter::parse_file(*in.xinfo);
        }
    }

    // Image source: CLI, then .xinfo, then XDS.INP.
    if (in.image) {
        m.image = {*in.image, MetadataSource::CommandLine};
    }
    for (const auto& [g, src] : {std::pair{&xinfo, MetadataSource::Xinfo},
                                 std::pair{&xds, MetadataSource::XdsInp}}) {
        if (!*g) {
            continue;
        }
        if (const auto p = image_from_import(**g)) {
            if (!m.image.set()) {
                m.image = {*p, src};
            } else if (std::filesystem::weakly_canonical(*p).parent_path() !=
                       std::filesystem::weakly_canonical(m.image.value).parent_path()) {
                m.warnings.push_back(
                    std::string("Image: ") + source_label(src, false) + " names " +
                    p->string() + " but " + m.image.value.string() + " (" +
                    source_label(m.image.source, false) + ") is used"
                );
            }
        }
    }
    if (!m.image.set()) {
        throw std::runtime_error("no image: give --image, or an --xinfo / --xds-inp naming one");
    }
    m.image.value = std::filesystem::absolute(m.image.value);
    m.hdf5 = is_container(m.image.value);
    if (!std::filesystem::exists(m.image.value)) {
        throw std::runtime_error("image not found: " + m.image.value.string());
    }

    // Frames on disk: an HDF5 master's count, or the consecutive CBF/SMV run.
    int first_on_disk = 1;
    if (m.hdf5) {
        m.frames_available = image::num_diffraction_images(m.image.value);
    } else {
        m.directory = m.image.value.parent_path();
        m.image_template = template_from_image_name(m.image.value.filename().string());
        first_on_disk = frame_number(m.image.value);
        int n = first_on_disk;
        std::error_code ec;
        while (n < first_on_disk + 1000000 &&
               std::filesystem::exists(frame_path(m.directory, m.image_template, n), ec)) {
            ++n;
        }
        m.frames_available = n - first_on_disk;
    }
    m.frame_range = {{first_on_disk, first_on_disk + std::max(m.frames_available, 1) - 1},
                     MetadataSource::ImageHeader};

    // Lowest priority: the first frame's header.
    const auto img = image::read_diffraction_image(m.image.value, m.hdf5 ? 1 : first_on_disk);
    const auto& md = img.metadata;
    if (md.wavelength_angstrom > 0.0) {
        m.wavelength_A = {md.wavelength_angstrom, MetadataSource::ImageHeader};
    }
    if (md.detector_distance_mm > 0.0) {
        m.distance_mm = {md.detector_distance_mm, MetadataSource::ImageHeader};
    }
    if (md.beam_x_mm > 0.0 || md.beam_y_mm > 0.0) {
        m.beam_mm = {{md.beam_x_mm, md.beam_y_mm}, MetadataSource::ImageHeader};
    }
    if (md.pixel_size_x_mm > 0.0) {
        m.pixel_mm = {{md.pixel_size_x_mm,
                       md.pixel_size_y_mm > 0.0 ? md.pixel_size_y_mm : md.pixel_size_x_mm},
                      MetadataSource::ImageHeader};
    }
    if (md.oscillation_range_deg > 0.0) {
        m.oscillation_deg = {md.oscillation_range_deg, MetadataSource::ImageHeader};
    }

    // Then XDS.INP, .xinfo and the command line, each per field.
    if (xds) {
        apply_import(m, *xds, MetadataSource::XdsInp);
    }
    if (xinfo) {
        apply_import(m, *xinfo, MetadataSource::Xinfo);
    }
    if (in.xparm) {
        const auto g = read_xparm(*in.xparm);
        apply_import(m, imported_from_xparm(g), in.xparm_source);
        if (g.has_amat) {
            std::array<double, 9> a{};
            for (int r = 0; r < 3; ++r) {
                for (int k = 0; k < 3; ++k) {
                    a[static_cast<std::size_t>(3 * r + k)] = g.amat(r, k);
                }
            }
            m.orientation = {a, in.xparm_source};
            m.orientation_control = xparm_orientation_lines(g);
        }
    }
    if (in.tilt_twist_mrad) {
        m.tilt_twist_mrad = {*in.tilt_twist_mrad, in.xparm_source};
    }
    if (in.range) {
        auto [a, b] = *in.range;
        if (b < a) {
            std::swap(a, b);
        }
        m.frame_range.assign({a, b}, MetadataSource::CommandLine);
    }
    const auto [f0, f1] = m.frame_range.value;
    const int last_on_disk = first_on_disk + m.frames_available - 1;
    if (f0 < first_on_disk || f1 > last_on_disk) {
        throw std::runtime_error(
            "frame range " + std::to_string(f0) + "–" + std::to_string(f1) + " (" +
            source_label(m.frame_range.source, m.hdf5) + ") is outside frames " +
            std::to_string(first_on_disk) + "–" + std::to_string(last_on_disk) +
            " of " + m.image.value.string()
        );
    }
    return m;
}

void write_xds_inp(const std::filesystem::path& path, const ExperimentModel& model) {
    std::ostringstream out;
    out << "! XDS.INP written by ProcessingLab from the resolved run options.\n";
    if (model.hdf5) {
        if (model.image.set()) {
            out << "NAME_TEMPLATE_OF_DATA_FRAMES=" << model.image.value.string() << '\n';
        }
    } else if (!model.image_template.empty()) {
        std::string tmpl = model.image_template;
        for (char& c : tmpl) {
            if (c == '#') {
                c = '?';
            }
        }
        const std::filesystem::path full =
            model.directory.empty() ? std::filesystem::path(tmpl) : model.directory / tmpl;
        out << "NAME_TEMPLATE_OF_DATA_FRAMES=" << full.string() << '\n';
    } else if (model.image.set()) {
        out << "NAME_TEMPLATE_OF_DATA_FRAMES=" << model.image.value.string() << '\n';
    }
    if (model.frame_range.set()) {
        out << "DATA_RANGE=" << model.frame_range.value.first << ' '
            << model.frame_range.value.second << '\n';
    }
    if (model.wavelength_A.set()) {
        out << "X-RAY_WAVELENGTH=" << fmt(model.wavelength_A.value, 6) << '\n';
    }
    if (model.distance_mm.set()) {
        out << "DETECTOR_DISTANCE=" << fmt(model.distance_mm.value, 3) << '\n';
    }
    if (model.pixel_mm.set()) {
        out << "QX=" << fmt(model.pixel_mm.value.first, 6)
            << " QY=" << fmt(model.pixel_mm.value.second, 6) << '\n';
    }
    if (model.beam_mm.set() && model.pixel_mm.set() &&
        model.pixel_mm.value.first > 0.0 && model.pixel_mm.value.second > 0.0) {
        const double orgx = model.beam_mm.value.first / model.pixel_mm.value.first;
        const double orgy = model.beam_mm.value.second / model.pixel_mm.value.second;
        out << "ORGX=" << fmt(orgx, 4) << " ORGY=" << fmt(orgy, 4) << '\n';
    }
    if (model.oscillation_deg.set()) {
        out << "OSCILLATION_RANGE=" << fmt(model.oscillation_deg.value, 4) << '\n';
    }
    if (model.cell.set()) {
        out << "UNIT_CELL_CONSTANTS=";
        for (double v : model.cell.value) {
            out << ' ' << fmt(v, 2);
        }
        out << '\n';
    }
    if (model.symmetry.set()) {
        out << "SPACE_GROUP_NUMBER=" << model.symmetry.value << '\n';
    }
    write_file_atomic(path, out.str());
}

std::vector<std::string> ExperimentModel::control_lines() const {
    std::vector<std::string> lines;
    if (!hdf5) {
        lines.push_back("DIRECTORY " + directory.string());
        lines.push_back("TEMPLATE " + image_template);
    }
    // Imported/CLI geometry only; header values are the session's defaults.
    ImportedSweepGeometry g = imported;
    g.has_directory = g.has_image = g.has_template = g.has_process_range = false;
    for (auto& l : g.to_control_lines()) {
        lines.push_back(std::move(l));
    }
    if (tilt_twist_mrad.set()) {
        lines.push_back("TILT " + fmt(tilt_twist_mrad.value.first, 4));
        lines.push_back("TWIST " + fmt(tilt_twist_mrad.value.second, 4));
    }
    lines.insert(lines.end(), orientation_control.begin(), orientation_control.end());
    return lines;
}

void ExperimentModel::write_summary(std::ostream& os) const {
    const auto tag = [&](MetadataSource s) {
        return std::string("  [") + source_label(s, hdf5) + "]";
    };
    os << "\n Experiment metadata\n";
    os << "  Image source:  " << image.value.string() << tag(image.source) << '\n';
    if (!hdf5) {
        os << "  Template:      " << (directory / image_template).string() << '\n';
    }
    if (wavelength_A.set()) {
        os << "  Wavelength:    " << fmt(wavelength_A.value, 5) << " Å" << tag(wavelength_A.source) << '\n';
    }
    if (distance_mm.set()) {
        os << "  Distance:      " << fmt(distance_mm.value, 2) << " mm" << tag(distance_mm.source) << '\n';
    }
    if (beam_mm.set()) {
        os << "  Beam:          " << fmt(beam_mm.value.first, 2) << " " << fmt(beam_mm.value.second, 2)
           << " mm" << tag(beam_mm.source) << '\n';
    }
    if (pixel_mm.set()) {
        os << "  Pixel:         " << fmt(pixel_mm.value.first, 4) << " " << fmt(pixel_mm.value.second, 4)
           << " mm" << tag(pixel_mm.source) << '\n';
    }
    if (oscillation_deg.set()) {
        os << "  Oscillation:   " << fmt(oscillation_deg.value, 4) << "°" << tag(oscillation_deg.source) << '\n';
    }
    if (cell.set()) {
        os << "  Cell:         ";
        for (double v : cell.value) {
            os << ' ' << fmt(v, 2);
        }
        os << tag(cell.source) << '\n';
    } else {
        os << "  Cell:          from autoindexing\n";
    }
    if (symmetry.set()) {
        os << "  Space group:   " << symmetry.value << tag(symmetry.source) << '\n';
    }
    if (tilt_twist_mrad.set()) {
        os << "  Tilt/twist:    " << fmt(tilt_twist_mrad.value.first, 2) << " " << fmt(tilt_twist_mrad.value.second, 2)
           << " mrad" << tag(tilt_twist_mrad.source) << '\n';
    }
    if (orientation.set()) {
        os << "  Orientation:   refined matrix, no autoindexing" << tag(orientation.source) << '\n';
    }
    os << "  Frame range:   " << frame_range.value.first << "–" << frame_range.value.second
       << tag(frame_range.source) << "  (" << frames_available << " on disk)\n";
    for (const auto& w : warnings) {
        os << "  *** WARNING: " << w << '\n';
    }
}

void apply_control_state(
    ExperimentModel& m, const control::ControlState& st, MetadataSource source
) {
    ImportedSweepGeometry g;
    const auto differs = [](double a, double b) { return rel_differs(a, b, 1e-9); };
    // Only values the script changed: imports already applied to ``st`` match.
    if (st.has_wavelength && !(m.wavelength_A.set() && !differs(m.wavelength_A.value, st.wavelength_A))) {
        g.has_wavelength = true;
        g.wavelength_A = st.wavelength_A;
    }
    if (st.has_distance && !(m.distance_mm.set() && !differs(m.distance_mm.value, st.distance_mm))) {
        g.has_distance = true;
        g.distance_mm = st.distance_mm;
    }
    if (st.has_beam && !(m.beam_mm.set() && !differs(m.beam_mm.value.first, st.beam_x_mm) &&
                         !differs(m.beam_mm.value.second, st.beam_y_mm))) {
        g.has_beam = true;
        g.beam_x_mm = st.beam_x_mm;
        g.beam_y_mm = st.beam_y_mm;
    }
    if (st.has_pixel && !(m.pixel_mm.set() && !differs(m.pixel_mm.value.first, st.pixel_x_mm) &&
                          !differs(m.pixel_mm.value.second, st.pixel_y_mm))) {
        g.has_pixel = true;
        g.pixel_x_mm = st.pixel_x_mm;
        g.pixel_y_mm = st.pixel_y_mm;
    }
    const std::array<double, 6> cell{st.cell.a, st.cell.b, st.cell.c,
                                     st.cell.alpha_deg, st.cell.beta_deg, st.cell.gamma_deg};
    if (st.has_cell && !(m.cell.set() && m.cell.value == cell)) {
        g.has_cell = true;
        g.cell_a = cell[0];
        g.cell_b = cell[1];
        g.cell_c = cell[2];
        g.cell_alpha = cell[3];
        g.cell_beta = cell[4];
        g.cell_gamma = cell[5];
    }
    if (st.has_symmetry && !(m.symmetry.set() && m.symmetry.value == st.symmetry_number)) {
        g.has_symmetry = true;
        g.symmetry_number = st.symmetry_number;
    }
    if (st.process_pending &&
        m.frame_range.value != std::pair<int, int>{st.process_first, st.process_last}) {
        g.has_process_range = true;
        g.process_first = st.process_first;
        g.process_last = st.process_last;
    }
    apply_import(m, g, source);
}

Json ExperimentModel::to_json() const {
    const auto field = [](const auto& f, auto&& encode, const char* unit) {
        Json j = Json::object();
        j["value"] = f.set() ? encode(f.value) : Json();
        if (unit != nullptr) {
            j["unit"] = unit;
        }
        j["source"] = source_id(f.source);
        Json ov = Json::array();
        for (const auto& [v, src] : f.overrides) {
            Json o = Json::object();
            o["value"] = encode(v);
            o["source"] = source_id(src);
            ov.push_back(std::move(o));
        }
        j["overrides"] = std::move(ov);
        return j;
    };
    const auto num = [](double v) { return Json(v); };
    const auto pair = [](const std::pair<double, double>& v) {
        return Json(Json::Array{Json(v.first), Json(v.second)});
    };
    Json e = Json::object();
    e["wavelength"] = field(wavelength_A, num, "angstrom");
    e["detector_distance"] = field(distance_mm, num, "mm");
    e["beam_centre"] = field(beam_mm, pair, "mm");
    e["pixel_size"] = field(pixel_mm, pair, "mm");
    e["oscillation_width"] = field(oscillation_deg, num, "deg");
    e["unit_cell"] = field(cell, [](const std::array<double, 6>& c) {
        Json::Array a;
        for (double v : c) {
            a.emplace_back(v);
        }
        return Json(std::move(a));
    }, "angstrom/deg");
    e["space_group_number"] = field(symmetry, [](int v) { return Json(v); }, nullptr);
    e["detector_tilt_twist"] = field(tilt_twist_mrad, pair, "mrad");
    e["orientation_matrix"] = field(orientation, [](const std::array<double, 9>& a) {
        Json::Array rows;
        for (std::size_t r = 0; r < 3; ++r) {
            rows.emplace_back(Json::Array{Json(a[3 * r]), Json(a[3 * r + 1]), Json(a[3 * r + 2])});
        }
        return Json(std::move(rows));
    }, nullptr);
    e["frame_range"] = field(frame_range, [](const std::pair<int, int>& r) {
        Json j = Json::object();
        j["first"] = r.first;
        j["last"] = r.second;
        return j;
    }, nullptr);
    return e;
}

std::filesystem::path create_output_directory(
    const std::optional<std::filesystem::path>& requested,
    const std::filesystem::path& parent,
    std::time_t now
) {
    namespace fs = std::filesystem;
    std::error_code ec;
    if (requested) {
        const fs::path dir = fs::absolute(*requested);
        if (dir.has_parent_path()) {
            fs::create_directories(dir.parent_path(), ec);
        }
        if (!fs::create_directory(dir, ec)) {
            throw std::runtime_error(
                ec ? "cannot create output directory " + dir.string() + ": " + ec.message()
                   : "output directory " + dir.string() +
                         " already exists; choose a new --output-dir (runs never reuse one)"
            );
        }
        return dir;
    }
    std::tm tm{};
    localtime_r(&now, &tm);
    char stamp[32];
    std::strftime(stamp, sizeof stamp, "processing_%Y%m%d_%H%M%S", &tm);
    for (int k = 1; k < 1000; ++k) {
        const fs::path dir = fs::absolute(parent) /
                             (k == 1 ? std::string(stamp) : std::string(stamp) + "_" + std::to_string(k));
        if (fs::create_directory(dir, ec)) {
            return dir;
        }
        if (ec) {
            throw std::runtime_error("cannot create output directory " + dir.string() + ": " + ec.message());
        }
    }
    throw std::runtime_error(std::string("no free output directory name for ") + stamp);
}

}  // namespace mosflm::io
