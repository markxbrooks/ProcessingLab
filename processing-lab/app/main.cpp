#include "mosflm/control/ControlSession.hpp"
#include "mosflm/crystal/MtzAdapter.hpp"
#include "mosflm/diffraction/ProcessPipeline.hpp"
#include "mosflm/diffraction/ProfileFitter.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/diffraction/SummationIntegrator.hpp"
#include "mosflm/geometry/UnitCell.hpp"
#include "mosflm/image/DiffractionImage.hpp"
#include "mosflm/legacy/LegacyMosflmAdapter.hpp"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace {

void print_usage(const char* argv0) {
    std::cerr
        << "Usage: " << argv0 << " <command> [args]\n"
        << "Commands:\n"
        << "  help                         Show this message\n"
        << "  version                      Print version / oracle hint\n"
        << "  predict-still [dump]         Analytical stills PREDICT → dump\n"
        << "  predict-osc [dump]           Analytical oscillation PREDICT → dump\n"
        << "  vertical-slice [mtz]         Predict→plant→integrate→MTZ smoke\n"
        << "  cbf-slice <cbf> [mtz]        Read miniCBF → predict → integrate → MTZ\n"
        << "  profile-fit [mtz]            Synthetic Gaussian peaks → profile fit → MTZ\n"
        << "  profile-fit-cbf <cbf> [mtz]  miniCBF → predict → profile fit → MTZ\n"
        << "  process-insulin [img] [mtz]  SMV insulin → predict → sum → MTZ (D1)\n"
        << "  run <script.in> [--image img] [mtz]  CONTROL keyword PROCESS (D5)\n";
}

int cmd_version() {
    std::cout << "processing-lab 0.1.0\n";
    std::cout << mosflm::legacy::LegacyMosflmAdapter::oracle_hint() << '\n';
    mosflm::UnitCell cell{58.0, 58.0, 150.0, 90.0, 90.0, 90.0};
    std::cout << "demo cell: " << cell.describe() << '\n';
    return 0;
}

int cmd_predict_still(const std::filesystem::path& dump) {
    const auto req =
        mosflm::diffraction::ReflectionPredictor::analytical_still_fixture();
    const auto result =
        mosflm::legacy::LegacyMosflmAdapter::predict_still(req);
    mosflm::diffraction::ReflectionPredictor::write_dump(
        dump, result.reflections
    );
    std::cout << "wrote " << result.reflections.size()
              << " stills predictions → " << dump << '\n';
    return 0;
}

int cmd_predict_osc(const std::filesystem::path& dump) {
    const auto req =
        mosflm::diffraction::ReflectionPredictor::analytical_oscillation_fixture();
    const auto result =
        mosflm::diffraction::ReflectionPredictor{}.predict(req);
    mosflm::diffraction::ReflectionPredictor::write_dump(
        dump, result.reflections
    );
    int n_partial = 0;
    for (const auto& s : result.reflections) {
        if (!s.full) {
            ++n_partial;
        }
    }
    std::cout << "wrote " << result.reflections.size() << " osc predictions ("
              << n_partial << " partial) → " << dump << '\n';
    return 0;
}

void plant_box(
    mosflm::image::DiffractionImage& image,
    int col,
    int row,
    int half,
    std::int32_t value
) {
    for (int dr = -half; dr <= half; ++dr) {
        for (int dc = -half; dc <= half; ++dc) {
            const int c = col + dc;
            const int r = row + dr;
            if (c < 0 || r < 0 || c >= image.metadata.width ||
                r >= image.metadata.height) {
                continue;
            }
            const std::size_t idx =
                static_cast<std::size_t>(r) *
                    static_cast<std::size_t>(image.metadata.width) +
                static_cast<std::size_t>(c);
            image.pixels[idx] += value;
        }
    }
}

void plant_gaussian(
    mosflm::image::DiffractionImage& image,
    int col0,
    int row0,
    double amplitude,
    double background,
    double sigma_px,
    int half
) {
    const double inv_2s2 = 1.0 / (2.0 * sigma_px * sigma_px);
    for (int dr = -half; dr <= half; ++dr) {
        for (int dc = -half; dc <= half; ++dc) {
            const int c = col0 + dc;
            const int r = row0 + dr;
            if (c < 0 || r < 0 || c >= image.metadata.width ||
                r >= image.metadata.height) {
                continue;
            }
            const double rr = static_cast<double>(dc * dc + dr * dr);
            const double w = std::exp(-rr * inv_2s2);
            const std::size_t idx =
                static_cast<std::size_t>(r) *
                    static_cast<std::size_t>(image.metadata.width) +
                static_cast<std::size_t>(c);
            image.pixels[idx] = static_cast<std::int32_t>(
                std::lround(background + amplitude * w)
            );
        }
    }
}

int cmd_vertical_slice(const std::filesystem::path& mtz_path) {
    auto req =
        mosflm::diffraction::ReflectionPredictor::analytical_still_fixture();
    req.settings.dmin_A = 8.0;
    const auto preds =
        mosflm::legacy::LegacyMosflmAdapter::predict_still(req).reflections;

    mosflm::image::DiffractionImage image;
    image.metadata.width = req.detector.nfast;
    image.metadata.height = req.detector.nslow;
    image.pixels.assign(
        static_cast<std::size_t>(image.metadata.width) *
            static_cast<std::size_t>(image.metadata.height),
        0
    );
    constexpr int half = 2;
    constexpr std::int32_t peak = 40;
    for (const auto& p : preds) {
        int col = 0;
        int row = 0;
        mosflm::diffraction::SummationIntegrator::mm_to_image_pixel(
            req.detector, p.x_mm, p.y_mm, col, row
        );
        plant_box(image, col, row, half, peak);
    }

    mosflm::diffraction::IntegrationSettings settings;
    settings.box_half_px = half;
    settings.bg_frame_px = 1;
    const mosflm::diffraction::IntegrationRequest ireq{
        image, req.detector, preds, settings};
    const auto integrated =
        mosflm::legacy::LegacyMosflmAdapter::integrate_summation(ireq);

    std::vector<mosflm::crystal::SimpleReflection> simple;
    for (const auto& r : integrated) {
        if (!r.ok) {
            continue;
        }
        simple.push_back(
            {r.h, r.k, r.l, static_cast<float>(r.intensity),
             static_cast<float>(r.sigma)}
        );
    }
    mosflm::crystal::MtzAdapter::write_simple(mtz_path, req.crystal, simple);
    std::cout << "vertical-slice: preds=" << preds.size()
              << " integrated=" << simple.size() << " → " << mtz_path << '\n';
    return 0;
}

int cmd_cbf_slice(
    const std::filesystem::path& cbf_path,
    const std::filesystem::path& mtz_path
) {
    const auto image = mosflm::image::CbfImageReader{}.read(cbf_path);
    auto req =
        mosflm::diffraction::ReflectionPredictor::analytical_still_fixture();
    req.detector = image.metadata.to_detector_geometry();
    req.detector.beam_x_mm = 0.0;
    req.detector.beam_y_mm = 0.0;
    req.scan.phi_start_deg = image.metadata.oscillation_start_deg;
    req.scan.phi_end_deg = image.metadata.oscillation_start_deg;
    req.settings.dmin_A = 4.0;

    const auto preds =
        mosflm::legacy::LegacyMosflmAdapter::predict_still(req).reflections;
    mosflm::diffraction::IntegrationSettings settings;
    settings.box_half_px = 2;
    settings.bg_frame_px = 1;
    const mosflm::diffraction::IntegrationRequest ireq{
        image, req.detector, preds, settings};
    const auto integrated =
        mosflm::legacy::LegacyMosflmAdapter::integrate_summation(ireq);

    std::vector<mosflm::crystal::SimpleReflection> simple;
    for (const auto& r : integrated) {
        if (!r.ok) {
            continue;
        }
        simple.push_back(
            {r.h, r.k, r.l, static_cast<float>(r.intensity),
             static_cast<float>(r.sigma)}
        );
    }
    mosflm::crystal::MtzAdapter::write_simple(mtz_path, req.crystal, simple);
    std::cout << "cbf-slice: " << cbf_path << " dims=" << image.metadata.width
              << "x" << image.metadata.height << " preds=" << preds.size()
              << " integrated=" << simple.size() << " → " << mtz_path << '\n';
    return 0;
}

int cmd_profile_fit(const std::filesystem::path& mtz_path) {
    auto req =
        mosflm::diffraction::ReflectionPredictor::analytical_still_fixture();
    req.settings.dmin_A = 8.0;
    const auto preds =
        mosflm::legacy::LegacyMosflmAdapter::predict_still(req).reflections;

    mosflm::image::DiffractionImage image;
    image.metadata.width = req.detector.nfast;
    image.metadata.height = req.detector.nslow;
    image.pixels.assign(
        static_cast<std::size_t>(image.metadata.width) *
            static_cast<std::size_t>(image.metadata.height),
        0
    );
    constexpr int half = 5;
    constexpr double amp = 80.0;
    constexpr double bg = 5.0;
    constexpr double sig = 1.5;
    for (const auto& p : preds) {
        int col = 0;
        int row = 0;
        mosflm::diffraction::SummationIntegrator::mm_to_image_pixel(
            req.detector, p.x_mm, p.y_mm, col, row
        );
        plant_gaussian(image, col, row, amp, bg, sig, half);
    }

    mosflm::diffraction::ProfileFitSettings settings;
    settings.box_half_px = half;
    settings.sigma_px = sig;
    settings.refine_sigma = true;
    settings.apply_lp = true;
    const mosflm::diffraction::ProfileFitRequest freq{
        image, req.detector, preds, settings};
    const auto fitted =
        mosflm::legacy::LegacyMosflmAdapter::integrate_profile(freq);

    std::vector<mosflm::crystal::SimpleReflection> simple;
    int n_ok = 0;
    for (const auto& r : fitted) {
        if (!r.ok) {
            continue;
        }
        ++n_ok;
        simple.push_back(
            {r.h, r.k, r.l, static_cast<float>(r.intensity),
             static_cast<float>(r.sigma)}
        );
    }
    mosflm::crystal::MtzAdapter::write_simple(mtz_path, req.crystal, simple);
    std::cout << "profile-fit: preds=" << preds.size() << " fitted=" << n_ok
              << " → " << mtz_path << '\n';
    return 0;
}

int cmd_profile_fit_cbf(
    const std::filesystem::path& cbf_path,
    const std::filesystem::path& mtz_path
) {
    const auto image = mosflm::image::CbfImageReader{}.read(cbf_path);
    auto req =
        mosflm::diffraction::ReflectionPredictor::analytical_still_fixture();
    req.detector = image.metadata.to_detector_geometry();
    req.detector.beam_x_mm = 0.0;
    req.detector.beam_y_mm = 0.0;
    req.scan.phi_start_deg = image.metadata.oscillation_start_deg;
    req.scan.phi_end_deg = image.metadata.oscillation_start_deg;
    req.settings.dmin_A = 4.0;

    const auto preds =
        mosflm::legacy::LegacyMosflmAdapter::predict_still(req).reflections;
    mosflm::diffraction::ProfileFitSettings settings;
    settings.box_half_px = 5;
    settings.sigma_px = 1.5;
    settings.refine_sigma = true;
    settings.apply_lp = true;
    const mosflm::diffraction::ProfileFitRequest freq{
        image, req.detector, preds, settings};
    const auto fitted =
        mosflm::legacy::LegacyMosflmAdapter::integrate_profile(freq);

    std::vector<mosflm::crystal::SimpleReflection> simple;
    for (const auto& r : fitted) {
        if (!r.ok) {
            continue;
        }
        simple.push_back(
            {r.h, r.k, r.l, static_cast<float>(r.intensity),
             static_cast<float>(r.sigma)}
        );
    }
    mosflm::crystal::MtzAdapter::write_simple(mtz_path, req.crystal, simple);
    std::cout << "profile-fit-cbf: " << cbf_path
              << " preds=" << preds.size() << " fitted=" << simple.size()
              << " → " << mtz_path << '\n';
    return 0;
}

std::filesystem::path default_insulin_img() {
    if (const char* env = std::getenv("MOSFLM_IMAGE")) {
        return env;
    }
    if (const char* ccp4 = std::getenv("CCP4")) {
        const auto p =
            std::filesystem::path(ccp4) / "examples/data/insulin_1_001.img";
        if (std::filesystem::exists(p)) {
            return p;
        }
    }
    return "/Applications/ccp4-9/examples/data/insulin_1_001.img";
}

int cmd_process_insulin(
    const std::filesystem::path& img_path,
    const std::filesystem::path& mtz_path
) {
    const auto image = mosflm::image::SmvImageReader{}.read(img_path);
    const auto req =
        mosflm::diffraction::ProcessPipeline::known_indexed_insulin(image);
    const auto result = mosflm::diffraction::ProcessPipeline{}.run(req);
    mosflm::diffraction::ProcessPipeline::write_mtz(
        mtz_path, req.crystal, result, req.detector.wavelength_A
    );
    std::cout << "process-insulin: " << img_path
              << " preds=" << result.predictions.size()
              << " integrated=" << result.mtz_reflections.size() << " → "
              << mtz_path << '\n';
    return 0;
}

int cmd_run_script(
    const std::filesystem::path& script,
    const std::optional<std::filesystem::path>& image,
    const std::filesystem::path& mtz_path
) {
    mosflm::control::ControlSession session;
    const auto runs = session.run_script(script, image);
    std::size_t n_proc = 0;
    for (const auto& r : runs) {
        if (!r.ran_process) {
            continue;
        }
        ++n_proc;
        mosflm::diffraction::ProcessPipeline::write_mtz(
            mtz_path, r.crystal, r.process, session.state().wavelength_A
        );
        std::cout << "run: " << script << " image=" << r.image_used
                  << " preds=" << r.process.predictions.size()
                  << " integrated=" << r.process.mtz_reflections.size()
                  << " → " << mtz_path << '\n';
    }
    if (n_proc == 0) {
        std::cerr << "run: no PROCESS GO executed (need MATRIX/CELL/PROCESS)\n";
        return 1;
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        // Default: brief version banner (backward compatible).
        return cmd_version();
    }
    const std::string_view cmd = argv[1];
    try {
        if (cmd == "help" || cmd == "-h" || cmd == "--help") {
            print_usage(argv[0]);
            return 0;
        }
        if (cmd == "version") {
            return cmd_version();
        }
        if (cmd == "predict-still") {
            const std::filesystem::path dump =
                (argc >= 3) ? argv[2] : "predict_still_out.dump";
            return cmd_predict_still(dump);
        }
        if (cmd == "predict-osc") {
            const std::filesystem::path dump =
                (argc >= 3) ? argv[2] : "predict_osc_out.dump";
            return cmd_predict_osc(dump);
        }
        if (cmd == "vertical-slice") {
            const std::filesystem::path mtz =
                (argc >= 3) ? argv[2] : "vertical_slice.mtz";
            return cmd_vertical_slice(mtz);
        }
        if (cmd == "cbf-slice") {
            if (argc < 3) {
                print_usage(argv[0]);
                return 2;
            }
            const std::filesystem::path mtz =
                (argc >= 4) ? argv[3] : "cbf_slice.mtz";
            return cmd_cbf_slice(argv[2], mtz);
        }
        if (cmd == "profile-fit") {
            const std::filesystem::path mtz =
                (argc >= 3) ? argv[2] : "profile_fit.mtz";
            return cmd_profile_fit(mtz);
        }
        if (cmd == "profile-fit-cbf") {
            if (argc < 3) {
                print_usage(argv[0]);
                return 2;
            }
            const std::filesystem::path mtz =
                (argc >= 4) ? argv[3] : "profile_fit_cbf.mtz";
            return cmd_profile_fit_cbf(argv[2], mtz);
        }
        if (cmd == "process-insulin") {
            const std::filesystem::path img =
                (argc >= 3) ? argv[2] : default_insulin_img();
            const std::filesystem::path mtz =
                (argc >= 4) ? argv[3] : "insulin_cpp_process.mtz";
            return cmd_process_insulin(img, mtz);
        }
        if (cmd == "run") {
            if (argc < 3) {
                print_usage(argv[0]);
                return 2;
            }
            const std::filesystem::path script = argv[2];
            std::optional<std::filesystem::path> image;
            std::filesystem::path mtz = "control_process.mtz";
            for (int i = 3; i < argc; ++i) {
                const std::string_view a = argv[i];
                if (a == "--image") {
                    if (i + 1 >= argc) {
                        std::cerr << "run: --image needs a path\n";
                        return 2;
                    }
                    image = argv[++i];
                } else {
                    mtz = argv[i];
                }
            }
            if (!image) {
                image = default_insulin_img();
            }
            return cmd_run_script(script, image, mtz);
        }
        std::cerr << "Unknown command: " << cmd << '\n';
        print_usage(argv[0]);
        return 2;
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << '\n';
        return 1;
    }
}
