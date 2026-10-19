#include "mosflm/legacy/MosflmCApi.h"
#include "test_util.hpp"

#include <vector>

int main() {
    mosflm_cell6 cell{};
    mosflm_detector_c det{};
    double phi = 0.0, dmin = 0.0, etol = 0.0;
    mosflm_analytical_still_fixture(&cell, &det, &phi, &dmin, &etol);

    // Count-only mode: out == NULL, cap == 0 → total count, SUCCESS.
    int n_total = -1;
    int rc = mosflm_predict_still(&cell, &det, nullptr, nullptr, phi, dmin, etol,
                                  nullptr, 0, &n_total);
    CHECK_EQ(rc, MOSFLM_SUCCESS);
    CHECK_EQ(n_total, 656);

    // Full buffer.
    std::vector<mosflm_predicted_c> out(static_cast<std::size_t>(n_total));
    int n = -1;
    rc = mosflm_predict_still(&cell, &det, nullptr, "P 1", phi, dmin, etol,
                              out.data(), n_total, &n);
    CHECK_EQ(rc, MOSFLM_SUCCESS);
    CHECK_EQ(n, n_total);

    // Truncation: small buffer → first cap written, BUFFER_TOO_SMALL.
    std::vector<mosflm_predicted_c> small(10);
    int n_small = -1;
    rc = mosflm_predict_still(&cell, &det, nullptr, nullptr, phi, dmin, etol,
                              small.data(), 10, &n_small);
    CHECK_EQ(rc, MOSFLM_ERROR_BUFFER_TOO_SMALL);
    CHECK_EQ(n_small, 10);
    // Truncated prefix matches the full run.
    CHECK_EQ(small[0].h, out[0].h);
    CHECK_EQ(small[0].k, out[0].k);
    CHECK_EQ(small[0].l, out[0].l);

    // Invalid argument: NULL cell → INVALID_ARGUMENT, out_n reset to 0.
    int n_bad = 12345;
    rc = mosflm_predict_still(nullptr, &det, nullptr, nullptr, phi, dmin, etol,
                              out.data(), n_total, &n_bad);
    CHECK_EQ(rc, MOSFLM_ERROR_INVALID_ARGUMENT);
    CHECK_EQ(n_bad, 0);

    // NULL detector is likewise rejected.
    int n_bad2 = 12345;
    rc = mosflm_predict_still(&cell, nullptr, nullptr, nullptr, phi, dmin, etol,
                              out.data(), n_total, &n_bad2);
    CHECK_EQ(rc, MOSFLM_ERROR_INVALID_ARGUMENT);
    CHECK_EQ(n_bad2, 0);

    // Oscillation fixture: 368 predictions, all partial for this fixture.
    mosflm_cell6 ocell{};
    mosflm_detector_c odet{};
    double p_start = 0.0, p_end = 0.0, odmin = 0.0, oetol = 0.0, pstep = 0.0;
    mosflm_analytical_osc_fixture(&ocell, &odet, &p_start, &p_end, &odmin,
                                  &oetol, &pstep);
    int n_osc = -1;
    rc = mosflm_predict_osc(&ocell, &odet, nullptr, nullptr, p_start, p_end,
                            odmin, oetol, pstep, nullptr, 0, &n_osc);
    CHECK_EQ(rc, MOSFLM_SUCCESS);
    CHECK_EQ(n_osc, 368);

    std::vector<mosflm_predicted_c> osc(static_cast<std::size_t>(n_osc));
    int n_osc2 = 0;
    rc = mosflm_predict_osc(&ocell, &odet, nullptr, nullptr, p_start, p_end,
                            odmin, oetol, pstep, osc.data(), n_osc, &n_osc2);
    CHECK_EQ(rc, MOSFLM_SUCCESS);
    CHECK_EQ(n_osc2, 368);
    int n_partial = 0;
    for (const auto& s : osc) {
        if (!s.full) {
            ++n_partial;
        }
    }
    CHECK_EQ(n_partial, 368);

    // Integrate: NULL detector rejected; out_n reset.
    int n_int = 7;
    rc = mosflm_integrate_summation(nullptr, nullptr, 0, 0, nullptr, 0, 2, 1, 0,
                                    nullptr, 0, &n_int);
    CHECK_EQ(rc, MOSFLM_ERROR_INVALID_ARGUMENT);
    CHECK_EQ(n_int, 0);

    return mosflm_test::finish("c_api_contract");
}
