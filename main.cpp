/**
 * Placeholder entry for IDE run configs.
 *
 * Default ProcessingLab binaries (Phase D / D6):
 *   cmake -S mosflm-cpp -B mosflm-cpp/build && cmake --build mosflm-cpp/build
 *   ./mosflm-cpp/build/ipmosflm          # thin CCP4-compatible CONTROL
 *   ./mosflm-cpp/build/mosflm_cpp_cli    # full CLI
 *
 * Fortran amalgamation (optional oracle):
 *   make fortran-oracle
 */
#include <iostream>

int main() {
    std::cout
        << "ProcessingLab: use bin/ipmosflm or mosflm_cpp_cli "
           "(no amalgamation by default).\n";
    return 0;
}
