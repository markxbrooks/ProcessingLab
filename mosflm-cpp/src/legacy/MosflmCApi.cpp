#include "mosflm/legacy/MosflmCApi.h"

#include "mosflm/crystal/CrystalMetadata.hpp"
#include "mosflm/diffraction/LpCorrection.hpp"
#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "mosflm/diffraction/SummationIntegrator.hpp"
#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/geometry/OrientationMatrix.hpp"
#include "mosflm/geometry/UnitCell.hpp"
#include "mosflm/image/DiffractionImage.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

namespace {

mosflm::UnitCell to_cell(const mosflm_cell6& c) {
    return mosflm::UnitCell{
        c.a, c.b, c.c, c.alpha_deg, c.beta_deg, c.gamma_deg
    };
}

mosflm::DetectorGeometry to_det(const mosflm_detector_c& d) {
    mosflm::DetectorGeometry g;
    g.distance_mm = d.distance_mm;
    g.beam_x_mm = d.beam_x_mm;
    g.beam_y_mm = d.beam_y_mm;
    g.pixel_size_x_mm = d.pixel_size_x_mm;
    g.pixel_size_y_mm = d.pixel_size_y_mm;
    g.nfast = d.nfast;
    g.nslow = d.nslow;
    g.wavelength_A = d.wavelength_A;
    return g;
}

}  // namespace

extern "C" void mosflm_analytical_still_fixture(
    mosflm_cell6* cell,
    mosflm_detector_c* det,
    double* phi_deg,
    double* dmin_A,
    double* ewald_tol
) {
    const auto req =
        mosflm::diffraction::ReflectionPredictor::analytical_still_fixture();
    if (cell) {
        cell->a = req.crystal.cell.a;
        cell->b = req.crystal.cell.b;
        cell->c = req.crystal.cell.c;
        cell->alpha_deg = req.crystal.cell.alpha;
        cell->beta_deg = req.crystal.cell.beta;
        cell->gamma_deg = req.crystal.cell.gamma;
    }
    if (det) {
        det->distance_mm = req.detector.distance_mm;
        det->beam_x_mm = req.detector.beam_x_mm;
        det->beam_y_mm = req.detector.beam_y_mm;
        det->pixel_size_x_mm = req.detector.pixel_size_x_mm;
        det->pixel_size_y_mm = req.detector.pixel_size_y_mm;
        det->nfast = req.detector.nfast;
        det->nslow = req.detector.nslow;
        det->wavelength_A = req.detector.wavelength_A;
    }
    if (phi_deg) {
        *phi_deg = req.phi_start_deg;
    }
    if (dmin_A) {
        *dmin_A = req.dmin_A;
    }
    if (ewald_tol) {
        *ewald_tol = req.ewald_tolerance_Ainv;
    }
}

extern "C" int mosflm_predict_still(
    const mosflm_cell6* cell,
    const mosflm_detector_c* det,
    double phi_deg,
    double dmin_A,
    double ewald_tol,
    mosflm_predicted_c* out,
    int out_cap,
    int* out_n
) {
    if (!cell || !det || !out || !out_n || out_cap <= 0) {
        return 1;
    }
    *out_n = 0;
    try {
        mosflm::diffraction::PredictRequest req;
        const mosflm::UnitCell uc = to_cell(*cell);
        req.crystal =
            mosflm::crystal::CrystalMetadata::from_mosflm_cell(uc, "P 1");
        req.orientation = mosflm::OrientationMatrix::from_cell_and_rotation(uc);
        req.detector = to_det(*det);
        req.phi_start_deg = phi_deg;
        req.phi_end_deg = phi_deg;
        req.dmin_A = dmin_A;
        req.ewald_tolerance_Ainv = ewald_tol > 0.0 ? ewald_tol : 0.0025;

        const auto spots =
            mosflm::diffraction::ReflectionPredictor{}.predict(req);
        const int n = static_cast<int>(
            std::min(spots.size(), static_cast<std::size_t>(out_cap))
        );
        for (int i = 0; i < n; ++i) {
            out[i].h = spots[static_cast<std::size_t>(i)].h;
            out[i].k = spots[static_cast<std::size_t>(i)].k;
            out[i].l = spots[static_cast<std::size_t>(i)].l;
            out[i].x_mm = spots[static_cast<std::size_t>(i)].x_mm;
            out[i].y_mm = spots[static_cast<std::size_t>(i)].y_mm;
            out[i].phi_begin_deg =
                spots[static_cast<std::size_t>(i)].phi_begin_deg;
            out[i].phi_end_deg = spots[static_cast<std::size_t>(i)].phi_end_deg;
            out[i].full = spots[static_cast<std::size_t>(i)].full ? 1 : 0;
        }
        *out_n = n;
        return (spots.size() > static_cast<std::size_t>(out_cap)) ? 2 : 0;
    } catch (...) {
        return 3;
    }
}

extern "C" double mosflm_lp_factor(
    const mosflm_detector_c* det,
    double x_mm,
    double y_mm
) {
    if (!det || !(det->distance_mm > 0.0)) {
        return 1.0;
    }
    const double dx = x_mm - det->beam_x_mm;
    const double dy = y_mm - det->beam_y_mm;
    const double r = std::sqrt(dx * dx + dy * dy);
    const double two_theta = std::atan2(r, det->distance_mm);
    const double c2 = std::cos(two_theta);
    const double s2 = std::sin(two_theta);
    const double P = 0.5 * (1.0 + c2 * c2);  // unpolarized
    const double L = (std::fabs(s2) > 1e-8) ? (1.0 / std::fabs(s2)) : 1.0;
    return L * P;
}

extern "C" int mosflm_integrate_summation(
    const mosflm_detector_c* det,
    const int* pixels,
    int nfast,
    int nslow,
    const mosflm_predicted_c* preds,
    int n_preds,
    int box_half_px,
    int bg_frame_px,
    int apply_lp,
    mosflm_integrated_c* out,
    int out_cap,
    int* out_n
) {
    if (!det || !pixels || !preds || !out || !out_n || nfast <= 0 || nslow <= 0 ||
        n_preds < 0 || out_cap <= 0) {
        return 1;
    }
    *out_n = 0;
    try {
        mosflm::image::DiffractionImage image;
        image.metadata.width = nfast;
        image.metadata.height = nslow;
        const std::size_t npix =
            static_cast<std::size_t>(nfast) * static_cast<std::size_t>(nslow);
        image.pixels.resize(npix);
        for (std::size_t i = 0; i < npix; ++i) {
            image.pixels[i] = static_cast<std::int32_t>(pixels[i]);
        }

        mosflm::DetectorGeometry g = to_det(*det);
        g.nfast = nfast;
        g.nslow = nslow;

        std::vector<mosflm::diffraction::PredictedReflection> pred_v;
        pred_v.reserve(static_cast<std::size_t>(n_preds));
        for (int i = 0; i < n_preds; ++i) {
            mosflm::diffraction::PredictedReflection p;
            p.h = preds[i].h;
            p.k = preds[i].k;
            p.l = preds[i].l;
            p.x_mm = preds[i].x_mm;
            p.y_mm = preds[i].y_mm;
            p.phi_begin_deg = preds[i].phi_begin_deg;
            p.phi_end_deg = preds[i].phi_end_deg;
            p.full = preds[i].full != 0;
            pred_v.push_back(p);
        }

        mosflm::diffraction::IntegrateRequest ireq;
        ireq.box_half_px = box_half_px;
        ireq.bg_frame_px = bg_frame_px;
        ireq.apply_lp = apply_lp != 0;

        const auto integrated =
            mosflm::diffraction::SummationIntegrator{}.integrate(
                image, g, pred_v, ireq
            );
        const int n = static_cast<int>(
            std::min(integrated.size(), static_cast<std::size_t>(out_cap))
        );
        for (int i = 0; i < n; ++i) {
            const auto& r = integrated[static_cast<std::size_t>(i)];
            out[i].h = r.h;
            out[i].k = r.k;
            out[i].l = r.l;
            out[i].x_mm = r.x_mm;
            out[i].y_mm = r.y_mm;
            out[i].intensity = r.intensity;
            out[i].sigma = r.sigma;
            out[i].ok = r.ok ? 1 : 0;
        }
        *out_n = n;
        return (integrated.size() > static_cast<std::size_t>(out_cap)) ? 2 : 0;
    } catch (...) {
        return 3;
    }
}
