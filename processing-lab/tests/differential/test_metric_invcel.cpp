#include "mosflm/geometry/Clcalc.hpp"
#include "mosflm/geometry/UnitCell.hpp"
#include "test_util.hpp"

#include <fstream>
#include <string>

using mosflm::UnitCell;

namespace {

void read_fixture(const char* path, UnitCell& cell, double& astar, double& bstar,
                  double& cstar, double& volume) {
    std::ifstream in(path);
    std::string key;
    while (in >> key) {
        if (key.empty() || key[0] == '#') {
            std::string rest;
            std::getline(in, rest);
            continue;
        }
        double value = 0.0;
        in >> value;
        if (key == "a") {
            cell.a = value;
        } else if (key == "b") {
            cell.b = value;
        } else if (key == "c") {
            cell.c = value;
        } else if (key == "alpha") {
            cell.alpha_deg = value;
        } else if (key == "beta") {
            cell.beta_deg = value;
        } else if (key == "gamma") {
            cell.gamma_deg = value;
        } else if (key == "astar") {
            astar = value;
        } else if (key == "bstar") {
            bstar = value;
        } else if (key == "cstar") {
            cstar = value;
        } else if (key == "volume") {
            volume = value;
        }
    }
}

}  // namespace

int main() {
    const char* env = std::getenv("MOSFLM_ANALYTICAL_FIXTURE");
    if (env == nullptr) {
        std::cerr << "MOSFLM_ANALYTICAL_FIXTURE not set\n";
        return 1;
    }

    UnitCell cell;
    double astar = 0.0, bstar = 0.0, cstar = 0.0, volume = 0.0;
    read_fixture(env, cell, astar, bstar, cstar, volume);

    CHECK(cell.is_valid());
    CHECK_CLOSE(cell.volume(), volume, 1e-3);

    // METRIC / INVCEL: reciprocal cell lengths against the fixture.
    const UnitCell recip = mosflm::cell_from_metric_tensor(
        mosflm::metric_tensor(cell)
    ).reciprocal();
    CHECK_CLOSE(recip.a, astar, 1e-12);
    CHECK_CLOSE(recip.b, bstar, 1e-12);
    CHECK_CLOSE(recip.c, cstar, 1e-12);

    // Direct roundtrip through CLCALC preserves the cell.
    const UnitCell rt =
        mosflm::cell_from_metric_tensor(mosflm::metric_tensor(cell));
    CHECK_CLOSE(rt.a, cell.a, 1e-9);
    CHECK_CLOSE(rt.c, cell.c, 1e-9);
    CHECK_CLOSE(rt.volume(), volume, 1e-3);

    return mosflm_test::finish("metric_invcel");
}
