/**
 * Placeholder entry for IDE run configs.
 *
 * Default ProcessingLab binaries (Phase D / E):
 *   cmake -S processing-lab -B processing-lab/build && cmake --build processing-lab/build
 *   ./processing-lab/build/ProcessingLab    # C++ CONTROL (ProcessingLab.exe on Windows)
 *   ./processing-lab/build/mosflm_cpp_cli    # full CLI
 *
 * Fortran amalgamation (optional oracle → bin/ipmosflm):
 *   make fortran-oracle
 */
#include <iostream>

int main() {
    std::cout
        << "ProcessingLab: use bin/ProcessingLab or mosflm_cpp_cli "
           "(F77 oracle: make fortran-oracle → bin/ipmosflm).\n";
    return 0;
}
