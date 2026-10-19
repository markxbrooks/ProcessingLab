/**
 * REEK E0: DIVH/DIVV widen the Ewald gate and increase predicted spot count.
 */

#include "mosflm/diffraction/ReflectionPredictor.hpp"
#include "test_util.hpp"

#include <iostream>

using mosflm::diffraction::ReflectionPredictor;

int main() {
    auto req = ReflectionPredictor::analytical_still_fixture();
    req.settings.mosaic_spread_deg = 0.0;
    req.settings.divh_deg = 0.0;
    req.settings.divv_deg = 0.0;
    req.settings.ewald_tolerance_Ainv = 0.001;

    const auto base = ReflectionPredictor{}.predict(req);
    req.settings.divh_deg = 0.5;
    req.settings.divv_deg = 0.5;
    const auto widened = ReflectionPredictor{}.predict(req);

    std::cout << "prediction_divergence: n0=" << base.reflections.size()
              << " n_div=" << widened.reflections.size() << '\n';

    CHECK(!base.reflections.empty());
    CHECK(widened.reflections.size() > base.reflections.size());

    // Ice-ring exclusion removes some spots when enabled.
    req.settings.divh_deg = 0.0;
    req.settings.divv_deg = 0.0;
    req.settings.exclude_ice_rings = false;
    const auto no_ice = ReflectionPredictor{}.predict(req);
    req.settings.exclude_ice_rings = true;
    const auto ice = ReflectionPredictor{}.predict(req);
    std::cout << "prediction_divergence: n_ice_off=" << no_ice.reflections.size()
              << " n_ice_on=" << ice.reflections.size() << '\n';
    CHECK(ice.reflections.size() <= no_ice.reflections.size());

    return mosflm_test::finish("prediction_divergence");
}
