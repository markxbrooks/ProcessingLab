cmake_minimum_required(VERSION 3.24)

project(MosflmCpp
    VERSION 0.1.0
    LANGUAGES C CXX
    DESCRIPTION "Progressive C++ reimplementation of MOSFLM computational core"
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)
include(FetchContent)

option(MOSFLM_WITH_CCTBX
    "Optional CCTBX validator (Python POC / future targets); OFF by default"
    OFF
)

FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE
)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING_SAVED "${BUILD_TESTING}")
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(eigen)
set(BUILD_TESTING "${BUILD_TESTING_SAVED}" CACHE BOOL "" FORCE)

# Gemmi: header-oriented crystallographic infrastructure (not MOSFLM algorithms).
FetchContent_Declare(
    gemmi
    GIT_REPOSITORY https://github.com/project-gemmi/gemmi.git
    GIT_TAG v0.6.7
    GIT_SHALLOW TRUE
)
FetchContent_GetProperties(gemmi)
if(NOT gemmi_POPULATED)
    FetchContent_Populate(gemmi)
endif()

find_package(ZLIB)
if(NOT ZLIB_FOUND)
    message(FATAL_ERROR "ZLIB required for Gemmi MTZ I/O")
endif()

# Vendored CBFlib (Option A) — sibling of mosflm-cpp under the ipmosflm tree.
set(MOSFLM_IPMOSFLM_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." CACHE PATH
    "ipmosflm repository root (CBF + DPS libs)")
set(MOSFLM_CBF_INCLUDE_DIR "${MOSFLM_IPMOSFLM_ROOT}/cbf/include")
set(MOSFLM_CBF_LIB "${MOSFLM_IPMOSFLM_ROOT}/cbf/lib/libcbf.a")
if(NOT EXISTS "${MOSFLM_CBF_LIB}")
    message(FATAL_ERROR "CBFlib not found at ${MOSFLM_CBF_LIB}")
endif()

# Prebuilt DPS Bravais (Cambridge 5-arg API).
set(MOSFLM_DPS_INDEX_LIB "${MOSFLM_IPMOSFLM_ROOT}/lib/dps_index.a")
set(MOSFLM_DPS_UTIL_LIB "${MOSFLM_IPMOSFLM_ROOT}/lib/libdps.a")
if(NOT EXISTS "${MOSFLM_DPS_INDEX_LIB}" OR NOT EXISTS "${MOSFLM_DPS_UTIL_LIB}")
    message(WARNING
        "DPS libraries missing (${MOSFLM_DPS_INDEX_LIB}); "
        "BravaisClassifier tests will be skipped until index/ is built")
    set(MOSFLM_HAS_DPS_BRAVAIS OFF)
else()
    set(MOSFLM_HAS_DPS_BRAVAIS ON)
endif()

add_library(mosflm_geometry
    src/geometry/UnitCell.cpp
    src/geometry/ReciprocalLattice.cpp
    src/geometry/OrientationMatrix.cpp
    src/geometry/DetectorGeometry.cpp
    src/legacy/LegacyMosflmAdapter.cpp
    src/image/CbfImageReader.cpp
    src/crystal/CrystalMetadata.cpp
    src/crystal/SpaceGroupInfo.cpp
    src/crystal/MtzAdapter.cpp
    src/diffraction/ReflectionPredictor.cpp
    ${gemmi_SOURCE_DIR}/src/mtz.cpp
    ${gemmi_SOURCE_DIR}/src/sprintf.cpp
    ${gemmi_SOURCE_DIR}/src/gz.cpp
)

target_include_directories(mosflm_geometry
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${gemmi_SOURCE_DIR}/include
    PRIVATE
        ${MOSFLM_CBF_INCLUDE_DIR}
)

target_link_libraries(mosflm_geometry
    PUBLIC
        Eigen3::Eigen
        ZLIB::ZLIB
    PRIVATE
        ${MOSFLM_CBF_LIB}
        m
)

target_compile_definitions(mosflm_geometry
    PUBLIC
        GEMMI_WRITE_IMPLEMENTATION
)

if(MOSFLM_HAS_DPS_BRAVAIS)
    add_library(mosflm_indexing
        src/indexing/BravaisClassifier.cpp
    )
    target_include_directories(mosflm_indexing
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/include
        PRIVATE
            ${MOSFLM_IPMOSFLM_ROOT}/util
            ${MOSFLM_IPMOSFLM_ROOT}/index
    )
    target_link_libraries(mosflm_indexing
        PUBLIC
            mosflm_geometry
        PRIVATE
            ${MOSFLM_DPS_INDEX_LIB}
            ${MOSFLM_DPS_UTIL_LIB}
            m
    )
endif()

add_executable(mosflm_cpp_cli
    app/main.cpp
)

target_link_libraries(mosflm_cpp_cli
    PRIVATE
        mosflm_geometry
)

if(MOSFLM_WITH_CCTBX)
    find_program(MOSFLM_CCTBX_PYTHON
        NAMES cctbx.python libtbx.python
        PATHS
            /Applications/phenix-1.20.1-4487/build/bin
            ENV PATH
    )
    if(MOSFLM_CCTBX_PYTHON)
        add_test(
            NAME cctbx_cell_poc
            COMMAND ${MOSFLM_CCTBX_PYTHON}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tools/cctbx_cell_poc.py
        )
    else()
        message(WARNING "MOSFLM_WITH_CCTBX=ON but cctbx.python not found")
    endif()
endif()

if(BUILD_TESTING)
    add_executable(test_unit_cell
        tests/unit/test_unit_cell.cpp
    )
    target_link_libraries(test_unit_cell PRIVATE mosflm_geometry)
    add_test(NAME unit_cell COMMAND test_unit_cell)

    add_executable(test_orientation
        tests/unit/test_orientation.cpp
    )
    target_link_libraries(test_orientation PRIVATE mosflm_geometry)
    add_test(NAME orientation COMMAND test_orientation)

    add_executable(test_detector
        tests/unit/test_detector.cpp
    )
    target_link_libraries(test_detector PRIVATE mosflm_geometry)
    add_test(NAME detector COMMAND test_detector)

    add_executable(test_metric_invcel
        tests/differential/test_metric_invcel.cpp
    )
    target_link_libraries(test_metric_invcel PRIVATE mosflm_geometry)
    add_test(NAME metric_invcel COMMAND test_metric_invcel)

    add_executable(test_gemmi_crystal
        tests/differential/test_gemmi_crystal.cpp
    )
    target_link_libraries(test_gemmi_crystal PRIVATE mosflm_geometry)
    add_test(NAME gemmi_crystal COMMAND test_gemmi_crystal)

    add_executable(test_cbf_image_reader
        tests/unit/test_cbf_image_reader.cpp
    )
    target_link_libraries(test_cbf_image_reader PRIVATE mosflm_geometry)
    add_test(
        NAME cbf_image_reader
        COMMAND test_cbf_image_reader
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    set_tests_properties(cbf_image_reader PROPERTIES
        ENVIRONMENT
            "MOSFLM_CBF_FIXTURE=${CMAKE_CURRENT_SOURCE_DIR}/testdata/cbf/sample_minicbf.cbf"
    )

    add_executable(test_predict_stub
        tests/differential/test_predict_stub.cpp
    )
    target_link_libraries(test_predict_stub PRIVATE mosflm_geometry)
    add_test(
        NAME predict_stub
        COMMAND test_predict_stub
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    if(MOSFLM_HAS_DPS_BRAVAIS)
        add_executable(test_bravais_classifier
            tests/differential/test_bravais_classifier.cpp
        )
        target_link_libraries(test_bravais_classifier PRIVATE mosflm_indexing)
        add_test(NAME bravais_classifier COMMAND test_bravais_classifier)
    endif()
endif()
