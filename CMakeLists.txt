cmake_minimum_required(VERSION 3.24)

# ProcessingLab — progressive C++ reimplementation of MOSFLM (legacy Autotools
# tree remains under mosflm/; computational port under mosflm-cpp/).
project(ProcessingLab
    VERSION 0.1.0
    LANGUAGES C CXX
    DESCRIPTION "ProcessingLab: MOSFLM computational core in C++"
)

option(PROCESSINGLAB_BUILD_CPP
    "Build mosflm-cpp library, CLI, and CTest suite"
    ON
)

# Tiny pointer binary for IDE run configs (real CLI: mosflm_cpp_cli).
add_executable(ProcessingLab main.cpp)

if(PROCESSINGLAB_BUILD_CPP)
    add_subdirectory(mosflm-cpp)
endif()
