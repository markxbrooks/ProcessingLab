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

option(MOSFLM_WITH_FORTRAN
    "Build Fortran iso_c_binding Phase C bridge + test (needs gfortran)"
    ON
)

# Phase C: Fortran → C ABI smoke (optional if no Fortran compiler).
set(MOSFLM_HAS_FORTRAN OFF)
if(MOSFLM_WITH_FORTRAN)
    include(CheckLanguage)
    check_language(Fortran)
    if(CMAKE_Fortran_COMPILER)
        enable_language(Fortran)
        set(MOSFLM_HAS_FORTRAN ON)
    else()
        message(WARNING
            "MOSFLM_WITH_FORTRAN=ON but no Fortran compiler found; "
            "skipping fortran/ bridge")
    endif()
endif()

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
    src/geometry/Rotmat.cpp
    src/geometry/Dstar.cpp
    src/geometry/Clcalc.cpp
    src/geometry/Matset.cpp
    src/legacy/LegacyMosflmAdapter.cpp
    src/legacy/MosflmCApiValidation.cpp
    src/legacy/MosflmCApi.cpp
    src/legacy/MosflmFortranAdapters.cpp
    src/image/CbfImageReader.cpp
    src/crystal/CrystalMetadata.cpp
    src/crystal/SpaceGroupInfo.cpp
    src/crystal/MtzAdapter.cpp
    src/diffraction/ReflectionPredictor.cpp
    src/diffraction/SummationIntegrator.cpp
    src/diffraction/ProfileFitter.cpp
    src/diffraction/LpCorrection.cpp
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
            ${MOSFLM_IPMOSFLM_ROOT}/include
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
    set_tests_properties(predict_stub PROPERTIES
        ENVIRONMENT
            "MOSFLM_PREDICT_ORACLE=${CMAKE_CURRENT_SOURCE_DIR}/testdata/oracle/predict_still.dump"
    )

    add_executable(test_summation_integrate
        tests/differential/test_summation_integrate.cpp
    )
    target_link_libraries(test_summation_integrate PRIVATE mosflm_geometry)
    add_test(
        NAME summation_integrate
        COMMAND test_summation_integrate
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_dstar_rotmat
        tests/unit/test_dstar_rotmat.cpp
    )
    target_link_libraries(test_dstar_rotmat PRIVATE mosflm_geometry)
    add_test(NAME dstar_rotmat COMMAND test_dstar_rotmat)

    add_executable(test_clcalc_matset
        tests/unit/test_clcalc_matset.cpp
    )
    target_link_libraries(test_clcalc_matset PRIVATE mosflm_geometry)
    add_test(NAME clcalc_matset COMMAND test_clcalc_matset)

    add_executable(test_vertical_slice
        tests/differential/test_vertical_slice.cpp
    )
    target_link_libraries(test_vertical_slice PRIVATE mosflm_geometry)
    add_test(NAME vertical_slice COMMAND test_vertical_slice)

    add_executable(test_predict_osc
        tests/differential/test_predict_osc.cpp
    )
    target_link_libraries(test_predict_osc PRIVATE mosflm_geometry)
    add_test(
        NAME predict_osc
        COMMAND test_predict_osc
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    set_tests_properties(predict_osc PROPERTIES
        ENVIRONMENT
            "MOSFLM_PREDICT_OSC_ORACLE=${CMAKE_CURRENT_SOURCE_DIR}/testdata/oracle/predict_osc.dump"
    )

    add_executable(test_cbf_vertical_slice
        tests/differential/test_cbf_vertical_slice.cpp
    )
    target_link_libraries(test_cbf_vertical_slice PRIVATE mosflm_geometry)
    add_test(
        NAME cbf_vertical_slice
        COMMAND test_cbf_vertical_slice
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    set_tests_properties(cbf_vertical_slice PROPERTIES
        ENVIRONMENT
            "MOSFLM_CBF_FIXTURE=${CMAKE_CURRENT_SOURCE_DIR}/testdata/cbf/sample_minicbf.cbf"
    )

    add_executable(test_profile_fitter
        tests/differential/test_profile_fitter.cpp
    )
    target_link_libraries(test_profile_fitter PRIVATE mosflm_geometry)
    add_test(NAME profile_fitter COMMAND test_profile_fitter)

    add_executable(test_c_api_lp
        tests/unit/test_c_api_lp.cpp
    )
    target_link_libraries(test_c_api_lp PRIVATE mosflm_geometry)
    add_test(NAME c_api_lp COMMAND test_c_api_lp)

    add_executable(test_c_api_contract
        tests/unit/test_c_api_contract.cpp
    )
    target_link_libraries(test_c_api_contract PRIVATE mosflm_geometry)
    add_test(NAME c_api_contract COMMAND test_c_api_contract)

    add_executable(test_ub_from_psi
        tests/unit/test_ub_from_psi.cpp
    )
    target_link_libraries(test_ub_from_psi PRIVATE mosflm_geometry)
    add_test(NAME ub_from_psi COMMAND test_ub_from_psi)

    if(MOSFLM_HAS_FORTRAN)
        add_executable(test_fortran_c_api
            fortran/mosflm_c_bind.f90
            fortran/test_fortran_c_api.f90
        )
        target_link_libraries(test_fortran_c_api PRIVATE mosflm_geometry)
        # C++ provides the ABI implementation; keep C++ as the link driver and
        # pull in the GNU Fortran runtime when needed.
        set_target_properties(test_fortran_c_api PROPERTIES
            LINKER_LANGUAGE CXX
            Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_mods
        )
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
            target_link_libraries(test_fortran_c_api PRIVATE gfortran)
        endif()
        add_test(NAME fortran_c_api COMMAND test_fortran_c_api)

        # Fixed-form F77 call into predict_still_cpp_ (amalgamation-style).
        add_executable(test_predict_still_f77
            fortran/test_predict_still_f77.f
        )
        set_source_files_properties(fortran/test_predict_still_f77.f
            PROPERTIES Fortran_FORMAT FIXED
        )
        target_link_libraries(test_predict_still_f77 PRIVATE mosflm_geometry)
        set_target_properties(test_predict_still_f77 PROPERTIES
            LINKER_LANGUAGE CXX
        )
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
            target_link_libraries(test_predict_still_f77 PRIVATE gfortran)
        endif()
        add_test(NAME predict_still_f77 COMMAND test_predict_still_f77)

        # Amalgamation-named PREDICT_CPP_STILL bridge.
        add_executable(test_predict_cpp_bridge
            fortran/predict_cpp_bridge.f90
            fortran/test_predict_cpp_bridge.f90
        )
        target_link_libraries(test_predict_cpp_bridge PRIVATE mosflm_geometry)
        set_target_properties(test_predict_cpp_bridge PROPERTIES
            LINKER_LANGUAGE CXX
        )
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
            target_link_libraries(test_predict_cpp_bridge PRIVATE gfortran)
        endif()
        add_test(NAME predict_cpp_bridge COMMAND test_predict_cpp_bridge)

        add_executable(test_predict_cpp_psi
            fortran/predict_cpp_bridge.f90
            fortran/test_predict_cpp_psi.f90
        )
        target_link_libraries(test_predict_cpp_psi PRIVATE mosflm_geometry)
        set_target_properties(test_predict_cpp_psi PROPERTIES
            LINKER_LANGUAGE CXX
        )
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
            target_link_libraries(test_predict_cpp_psi PRIVATE gfortran)
        endif()
        add_test(NAME predict_cpp_psi COMMAND test_predict_cpp_psi)
    endif()

    add_executable(test_known_indexed_oracle
        tests/differential/test_known_indexed_oracle.cpp
    )
    target_link_libraries(test_known_indexed_oracle PRIVATE mosflm_geometry)
    add_test(
        NAME known_indexed_oracle
        COMMAND test_known_indexed_oracle
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
