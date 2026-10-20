cmake_minimum_required(VERSION 3.24)

# ProcessingLab — progressive C++ reimplementation of MOSFLM (legacy Autotools
# tree remains under mosflm/; computational port under processing-lab/).
project(ProcessingLab
    VERSION 0.1.0
    LANGUAGES C CXX
    DESCRIPTION "ProcessingLab: MOSFLM computational core in C++"
)

option(PROCESSINGLAB_BUILD_CPP
    "Build processing-lab library, CLI, and CTest suite"
    ON
)

if(PROCESSINGLAB_BUILD_CPP)
    # Real binary: processing-lab/app/ipmosflm_main.cpp → target ProcessingLab.
    add_subdirectory(processing-lab)
else()
    # Tiny pointer binary for IDE run configs when the C++ port is disabled.
    add_executable(ProcessingLab main.cpp)
endif()
