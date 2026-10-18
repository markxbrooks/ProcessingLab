/**
 * Placeholder entry for IDE run configs.
 *
 * Default ProcessingLab binaries (Phase D / E):
 *   cmake -S mosflm-cpp -B mosflm-cpp/build && cmake --build mosflm-cpp/build
 *   ./mosflm-cpp/build/IntegrationLab    # C++ CONTROL (IntegrationLab.exe on Windows)
 *   ./mosflm-cpp/build/mosflm_cpp_cli    # full CLI
 *
 * Fortran amalgamation (optional oracle → bin/ipmosflm):
 *   make fortran-oracle
 */
#include <iostream>

int main() {
    std::cout
        << "ProcessingLab: use bin/IntegrationLab or mosflm_cpp_cli "
           "(F77 oracle: make fortran-oracle → bin/ipmosflm).\n";
    return 0;
}
