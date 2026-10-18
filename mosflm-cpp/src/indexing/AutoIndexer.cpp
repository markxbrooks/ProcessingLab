#include "mosflm/indexing/AutoIndexer.hpp"

#include "mosflm/indexing/BravaisSelect.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

extern "C" {
#include "util.h"
#include "index_def.h"
#include "cryst_matrix.h"
}

namespace mosflm::indexing {
namespace {

extern "C" void index_fortran_(
    int* numvect,
    int* numvectref,
    int* maxcyc,
    int* inicyc,
    float* convrg,
    int* numormat,
    float* stepscal,
    int* numrlp,
    int* numsptdps,
    int* npkfftmin,
    int* nmultfft,
    int* idumpth,
    float* theta1d,
    float* theta2d,
    float* sinthmin,
    float* volmin,
    float* hkldev,
    float* hklmin,
    int dumpfft[],
    int* dumpvec,
    int* dumpref,
    int* dumpgen,
    int* dumpch,
    int* number_xyz,
    int* nlattice,
    float* x,
    float* y,
    float* z,
    int* usespot,
    float* d_max,
    float* d_min,
    int* index_mode,
    float* a,
    float* wave,
    int* ierrflg,
    float* tarcell,
    int* centering
);

UnitCell cell_from_dps(const struct cell_struct& c) {
    UnitCell out;
    out.a = c.a;
    out.b = c.b;
    out.c = c.c;
    constexpr double kRad = 180.0 / 3.14159265358979323846;
    out.alpha_deg = c.alpha * kRad;
    out.beta_deg = c.beta * kRad;
    out.gamma_deg = c.gamma * kRad;
    return out;
}

Eigen::Matrix3d amat_from_dps_om(const struct matrix_3x3& om, double wave) {
    // ind_fortran.c: MOSFLM A = λ · permutation of DPS OM rows.
    Eigen::Matrix3d A;
    A(0, 0) = wave * om.r[0].z;
    A(0, 1) = wave * om.r[0].x;
    A(0, 2) = wave * om.r[0].y;
    A(1, 0) = wave * om.r[1].z;
    A(1, 1) = wave * om.r[1].x;
    A(1, 2) = wave * om.r[1].y;
    A(2, 0) = wave * om.r[2].z;
    A(2, 1) = wave * om.r[2].x;
    A(2, 2) = wave * om.r[2].y;
    return A;
}

struct matrix_3x3 eigen_to_dps(const Eigen::Matrix3d& M) {
    struct matrix_3x3 om{};
    for (int i = 0; i < 3; ++i) {
        om.r[i].x = static_cast<float>(M(i, 0));
        om.r[i].y = static_cast<float>(M(i, 1));
        om.r[i].z = static_cast<float>(M(i, 2));
    }
    return om;
}

}  // namespace

AutoIndexResult AutoIndexer::index(
    const SpotList& spots,
    const AutoIndexGeometry& geometry,
    const AutoIndexSettings& settings
) const {
    AutoIndexResult result;
    if (spots.spots.empty()) {
        return result;
    }

    const double beam_x =
        (geometry.beam_x_mm > 0.0) ? geometry.beam_x_mm : spots.header.beam_x_mm;
    const double beam_y =
        (geometry.beam_y_mm > 0.0) ? geometry.beam_y_mm : spots.header.beam_y_mm;

    std::vector<std::size_t> keep;
    keep.reserve(spots.spots.size());
    for (std::size_t i = 0; i < spots.spots.size(); ++i) {
        const auto& s = spots.spots[i];
        const double sig = (s.sigma > 1e-6) ? s.sigma : 1.0;
        if (s.intensity / sig >= static_cast<double>(settings.min_isigi)) {
            keep.push_back(i);
        }
    }
    const int n = static_cast<int>(keep.size());
    if (n > N_RS) {
        throw std::runtime_error("AutoIndexer: too many spots for DPS N_RS");
    }
    if (n < 20) {
        result.ierr = -2;
        result.n_spots_used = static_cast<std::size_t>(n);
        return result;
    }

    // MOSFLM DPSCHOOSE: FR_MAT (OMEGA + INVERTX) → triangular xyz → ROTMATY.
    // .spt INVFLAG≠0 ⇒ FINVERTX=-1 (ADSC). ROTMATY takes φ in degrees (DTOR).
    // Spot-line field 4 is φ_end (insulin: 0.5°). At TWOTHETA=MU=0 the γ/ν
    // path equals this triangular form.
    constexpr double kDeg2Rad = 3.14159265358979323846 / 180.0;
    const double dist = geometry.distance_mm;
    const double inv_l = 1.0 / geometry.wavelength_A;
    const double omega = spots.header.omega_deg * kDeg2Rad;
    const double finvertx = (spots.header.invflag != 0) ? -1.0 : 1.0;
    const double fr11 = std::cos(omega) * finvertx;
    const double fr12 = std::sin(omega);
    const double fr21 = -std::sin(omega) * finvertx;
    const double fr22 = std::cos(omega);
    const int n_xyz = n;
    result.n_spots_used = static_cast<std::size_t>(n_xyz);

    std::vector<float> x(static_cast<std::size_t>(n_xyz));
    std::vector<float> y(static_cast<std::size_t>(n_xyz));
    std::vector<float> z(static_cast<std::size_t>(n_xyz));
    std::vector<int> usespot(static_cast<std::size_t>(n_xyz), 1);
    for (int i = 0; i < n_xyz; ++i) {
        const auto& s = spots.spots[keep[static_cast<std::size_t>(i)]];
        const double xt = s.x_mm - beam_x;
        const double yt = s.y_mm - beam_y;
        const double xs = xt * fr11 + yt * fr12;
        const double ys = xt * fr21 + yt * fr22;
        const double rlen = std::sqrt(xs * xs + ys * ys + dist * dist);
        const double f = inv_l / rlen;
        const double qx = f * xs;
        const double qy = f * ys;
        const double qz = f * dist - inv_l;
        const double phi_rad = s.phi_end_deg * kDeg2Rad;
        const double c = std::cos(phi_rad);
        const double sn = std::sin(phi_rad);
        // ROTMATY (TRANSP∘RTMATS about Y): x'=c x − s z, z'=s x + c z
        x[static_cast<std::size_t>(i)] = static_cast<float>(c * qx - sn * qz);
        y[static_cast<std::size_t>(i)] = static_cast<float>(qy);
        z[static_cast<std::size_t>(i)] = static_cast<float>(sn * qx + c * qz);
    }

    int numvect = settings.nvect;
    int numvectref = settings.nvref;
    int maxcyc = settings.maxcyc;
    int inicyc = settings.inicyc;
    float convrg = settings.convrg;
    int numormat = settings.numormat;
    float stepscal = settings.stepscal;
    int numrlp = settings.numrlp;
    int numsptdps = settings.numsptdps;
    int npkfftmin = settings.npkfftmin;
    int nmultfft = settings.nmultfft;
    int idumpth = 0;
    float theta1d = settings.theta1d;
    float theta2d = settings.theta2d;
    float sinthmin = settings.sinthmin;
    float volmin = settings.volmin;
    float hkldev = settings.hkldev;
    float hklmin = settings.hklmin;
    int dumpfft[3] = {0, 0, 0};
    int dumpvec = 0;
    int dumpref = 0;
    int dumpgen = 0;
    int dumpch = 0;
    int number_xyz = n_xyz;
    int nlattice = 1;
    float d_max = static_cast<float>(geometry.dmax_A * settings.dmaxscal);
    float d_min = 0.0f;  // UMINCELL
    // MOSFLM DPSIND passes literal 1 (INDEX_AUTOMATIC is 0; unused in do_index).
    int index_mode = 1;
    float a9[9] = {};
    float wave = static_cast<float>(geometry.wavelength_A);
    int ierr = 0;
    float tarcell[6] = {0, 0, 0, 0, 0, 0};
    int centering = 0;

    index_fortran_(
        &numvect,
        &numvectref,
        &maxcyc,
        &inicyc,
        &convrg,
        &numormat,
        &stepscal,
        &numrlp,
        &numsptdps,
        &npkfftmin,
        &nmultfft,
        &idumpth,
        &theta1d,
        &theta2d,
        &sinthmin,
        &volmin,
        &hkldev,
        &hklmin,
        dumpfft,
        &dumpvec,
        &dumpref,
        &dumpgen,
        &dumpch,
        &number_xyz,
        &nlattice,
        x.data(),
        y.data(),
        z.data(),
        usespot.data(),
        &d_max,
        &d_min,
        &index_mode,
        a9,
        &wave,
        &ierr,
        tarcell,
        &centering
    );

    result.ierr = ierr;
    if (ierr != 0) {
        return result;
    }

    result.amat << a9[0], a9[1], a9[2], a9[3], a9[4], a9[5], a9[6], a9[7],
        a9[8];

    // Recover DPS OM from AMAT packing for cell + Bravais.
    struct matrix_3x3 om{};
    const double inv_w = 1.0 / geometry.wavelength_A;
    om.r[0].z = static_cast<float>(a9[0] * inv_w);
    om.r[0].x = static_cast<float>(a9[1] * inv_w);
    om.r[0].y = static_cast<float>(a9[2] * inv_w);
    om.r[1].z = static_cast<float>(a9[3] * inv_w);
    om.r[1].x = static_cast<float>(a9[4] * inv_w);
    om.r[1].y = static_cast<float>(a9[5] * inv_w);
    om.r[2].z = static_cast<float>(a9[6] * inv_w);
    om.r[2].x = static_cast<float>(a9[7] * inv_w);
    om.r[2].y = static_cast<float>(a9[8] * inv_w);

    struct cell_struct cell{};
    dircell_from_om(om, &cell);
    result.cell = cell_from_dps(cell);

    Eigen::Matrix3d om_e;
    for (int i = 0; i < 3; ++i) {
        om_e(i, 0) = om.r[i].x;
        om_e(i, 1) = om.r[i].y;
        om_e(i, 2) = om.r[i].z;
    }
    result.bravais = BravaisClassifier{}.classify(om_e, 1);

    const std::size_t pick = select_suggested_bravais(
        result.bravais,
        settings.max_penalty,
        settings.prefer_character,
        settings.prefer_system
    );
    if (pick != static_cast<std::size_t>(-1) && pick < result.bravais.size()) {
        const auto& c = result.bravais[pick];
        result.selected_character = c.character;
        result.selected_system = c.system_id;
        result.selected_symmetry_number =
            symmetry_number_for_lattice(c.system_id);
        // Always apply Bravais OM → CELL / AMAT (SUGGLATT + prefer_*).
        const struct matrix_3x3 om_pick = eigen_to_dps(c.orientation);
        result.amat = amat_from_dps_om(om_pick, geometry.wavelength_A);
        struct cell_struct cell2{};
        dircell_from_om(om_pick, &cell2);
        result.cell = cell_from_dps(cell2);
    }

    result.orientation =
        OrientationMatrix::from_mosflm_orientation_matrix(result.amat);
    result.ok = result.cell.is_valid();
    return result;
}

}  // namespace mosflm::indexing
