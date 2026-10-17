/**
 * Placeholder entry — the ProcessingLab CLI lives in mosflm-cpp
 * (``mosflm_cpp_cli``). Prefer:
 *
 *   cmake -S mosflm-cpp -B mosflm-cpp/build && cmake --build mosflm-cpp/build
 *   ./mosflm-cpp/build/mosflm_cpp_cli version
 *
 * Legacy Autotools binary remains ``bin/ipmosflm`` (MOSFLM keyword interface).
 */
#include <iostream>

int main() {
    std::cout
        << "ProcessingLab: use mosflm_cpp_cli (mosflm-cpp/) or bin/ipmosflm.\n";
    return 0;
}
