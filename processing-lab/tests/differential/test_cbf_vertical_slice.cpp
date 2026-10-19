#include "mosflm/crystal/MtzAdapter.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/diffraction/SummationIntegrator.hpp"
#include "mosflm/image/DiffractionImage.hpp"
#include "mosflm/legacy/LegacyMosflmAdapter.hpp"
#include "test_util.hpp"

#include <filesystem>
#include <vector>

using mosflm::crystal::MtzAdapter;
using mosflm::crystal::SimpleReflection;
using mosflm::diffraction::IntegrationRequest;
using mosflm::diffraction::IntegrationSettings;
using mosflm::diffraction::ReflectionPredictor;
using mosflm::diffraction::SummationIntegrator;
using mosflm::image::CbfImageReader;

int main() {
    const char* env = std::getenv("MOSFLM_CBF_FIXTURE");
    if (env == nullptr) {
        std::cerr << "MOSFLM_CBF_FIXTURE not set\n";
        return 1;
    }
    const auto image = CbfImageReader{}.read(std::filesystem::path(env));
    CHECK_EQ(image.metadata.width, 128);
    CHECK_EQ(image.metadata.height, 128);

    auto req = ReflectionPredictor::analytical_still_fixture();
    req.detector = image.metadata.to_detector_geometry();
    req.scan.phi_start_deg = image.metadata.oscillation_start_deg;
    req.scan.phi_end_deg = image.metadata.oscillation_start_deg;
    req.settings.dmin_A = 10.0;

    const auto preds =
        mosflm::legacy::LegacyMosflmAdapter::predict_still(req).reflections;
    CHECK(!preds.empty());

    IntegrationSettings settings;
    settings.box_half_px = 2;
    settings.bg_frame_px = 1;
    const IntegrationRequest ireq{image, req.detector, preds, settings};
    const auto integrated =
        mosflm::legacy::LegacyMosflmAdapter::integrate_summation(ireq);
    CHECK_EQ(integrated.size(), preds.size());

    std::vector<SimpleReflection> simple;
    for (const auto& r : integrated) {
        if (!r.ok) {
            continue;
        }
        simple.push_back({r.h, r.k, r.l, static_cast<float>(r.intensity),
                          static_cast<float>(r.sigma)});
    }

    const std::filesystem::path mtz =
        std::filesystem::temp_directory_path() / "mosflm_cbf_slice.mtz";
    MtzAdapter::write_simple(mtz, req.crystal, simple);
    CHECK(std::filesystem::exists(mtz));
    const auto back = MtzAdapter::read_simple(mtz);
    CHECK_EQ(back.size(), simple.size());

    std::filesystem::remove(mtz);
    return mosflm_test::finish("cbf_vertical_slice");
}
