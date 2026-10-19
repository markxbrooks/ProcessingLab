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

option(MOSFLM_FORTRAN_ORACLE
    "Build mosflm_predict_hook for Autotools amalgamation splice (opt-in; D6)"
    OFF
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

# Only add Fortran test targets when their driver sources are present.
set(MOSFLM_FORTRAN_TESTS OFF)
if(MOSFLM_HAS_FORTRAN AND
   EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/fortran/test_fortran_c_api.f90" AND
   EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/fortran/test_predict_still_f77.f" AND
   EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/fortran/test_predict_cpp_bridge.f90")
    set(MOSFLM_FORTRAN_TESTS ON)
else()
    message(WARNING
        "Fortran test drivers missing under fortran/; "
        "skipping Fortran test targets (Phase C hook library still built)")
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

# Vendored CBFlib (Option A) — sibling of processing-lab under the ipmosflm tree.
set(MOSFLM_IPMOSFLM_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." CACHE PATH
    "ProcessingLab repository root (CBF + DPS libs)")
set(MOSFLM_CBF_INCLUDE_DIR "${MOSFLM_IPMOSFLM_ROOT}/cbf/include")
set(MOSFLM_CBF_LIB "${MOSFLM_IPMOSFLM_ROOT}/cbf/lib/libcbf.a")
if(NOT EXISTS "${MOSFLM_CBF_LIB}")
    message(FATAL_ERROR "CBFlib not found at ${MOSFLM_CBF_LIB}")
endif()

# Prebuilt DPS (Cambridge Bravais + FFT indexer).
set(MOSFLM_DPS_INDEX_LIB "${MOSFLM_IPMOSFLM_ROOT}/lib/dps_index.a")
set(MOSFLM_DPS_UTIL_LIB "${MOSFLM_IPMOSFLM_ROOT}/lib/libdps.a")
set(MOSFLM_DPS_INDEXER_LIB "${MOSFLM_IPMOSFLM_ROOT}/lib/libdps_index.a")
set(MOSFLM_DPS_FFTPACK_LIB "${MOSFLM_IPMOSFLM_ROOT}/lib/libfftpack.a")
if(NOT EXISTS "${MOSFLM_DPS_INDEX_LIB}" OR NOT EXISTS "${MOSFLM_DPS_UTIL_LIB}")
    message(WARNING
        "DPS libraries missing (${MOSFLM_DPS_INDEX_LIB}); "
        "BravaisClassifier / AutoIndexer tests skipped until index/ is built")
    set(MOSFLM_HAS_DPS_BRAVAIS OFF)
else()
    set(MOSFLM_HAS_DPS_BRAVAIS ON)
endif()
set(MOSFLM_HAS_DPS_INDEXER OFF)
if(MOSFLM_HAS_DPS_BRAVAIS AND EXISTS "${MOSFLM_DPS_INDEXER_LIB}")
    set(MOSFLM_HAS_DPS_INDEXER ON)
endif()

# Optional HDF5 for Eiger NeXus / *_master.h5 (Phase E / E5 decode).
find_package(HDF5 COMPONENTS C HL)
set(MOSFLM_HAS_HDF5 OFF)
if(HDF5_FOUND)
    set(MOSFLM_HAS_HDF5 ON)
    message(STATUS "HDF5 found — enabling Hdf5ImageReader")
else()
    message(WARNING "HDF5 not found — Hdf5ImageReader stub only")
endif()

set(MOSFLM_GEOMETRY_SOURCES
    src/geometry/UnitCell.cpp
    src/geometry/CellFix.cpp
    src/geometry/ReciprocalLattice.cpp
    src/geometry/OrientationMatrix.cpp
    src/geometry/DetectorGeometry.cpp
    src/geometry/LabFrame.cpp
    src/geometry/Rotmat.cpp
    src/geometry/Dstar.cpp
    src/geometry/Clcalc.cpp
    src/geometry/Matset.cpp
    src/legacy/LegacyMosflmAdapter.cpp
    src/legacy/MosflmCApiValidation.cpp
    src/legacy/MosflmCApi.cpp
    src/legacy/MosflmFortranAdapters.cpp
    src/image/CbfImageReader.cpp
    src/image/SmvImageReader.cpp
    src/image/Hdf5ImageReader.cpp
    src/image/ImageIO.cpp
    src/imosflm/ImosflmSocket.cpp
    src/crystal/CrystalMetadata.cpp
    src/crystal/SpaceGroupInfo.cpp
    src/crystal/MtzAdapter.cpp
    src/diffraction/ReflectionPredictor.cpp
    src/diffraction/SummationIntegrator.cpp
    src/diffraction/ProcessPipeline.cpp
    src/diffraction/MeasurementBox.cpp
    src/diffraction/TrustedRegion.cpp
    src/diffraction/DetectorInit.cpp
    src/diffraction/StandardProfile.cpp
    src/diffraction/ProfileFitter.cpp
    src/diffraction/LpCorrection.cpp
    src/diffraction/DataQualityStats.cpp
    src/refinement/PostRefinement.cpp
    src/control/ControlSession.cpp
    src/io/ImportedSweepGeometry.cpp
    src/io/XinfoImporter.cpp
    src/io/XdsInpImporter.cpp
    src/io/XparmIO.cpp
    ${gemmi_SOURCE_DIR}/src/mtz.cpp
    ${gemmi_SOURCE_DIR}/src/sprintf.cpp
    ${gemmi_SOURCE_DIR}/src/gz.cpp
)

if(MOSFLM_HAS_HDF5)
    list(APPEND MOSFLM_GEOMETRY_SOURCES
        ${MOSFLM_IPMOSFLM_ROOT}/lz4/lz4.c
        ${MOSFLM_IPMOSFLM_ROOT}/lz4/h5zlz4.c
        ${MOSFLM_IPMOSFLM_ROOT}/bitshuffle/bitshuffle.c
        ${MOSFLM_IPMOSFLM_ROOT}/bitshuffle/bshuf_h5filter.c
    )
endif()

add_library(mosflm_geometry ${MOSFLM_GEOMETRY_SOURCES})

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

if(MOSFLM_HAS_HDF5)
    target_include_directories(mosflm_geometry PRIVATE
        ${HDF5_INCLUDE_DIRS}
        ${MOSFLM_IPMOSFLM_ROOT}/lz4
        ${MOSFLM_IPMOSFLM_ROOT}/bitshuffle
    )
    # ${HDF5_LIBRARIES} is often only libhdf5; H5LT lives in libhdf5_hl.
    target_link_libraries(mosflm_geometry PRIVATE ${HDF5_LIBRARIES})
    if(TARGET hdf5::hdf5_hl)
        target_link_libraries(mosflm_geometry PRIVATE hdf5::hdf5_hl)
    elseif(HDF5_HL_LIBRARIES)
        target_link_libraries(mosflm_geometry PRIVATE ${HDF5_HL_LIBRARIES})
    else()
        target_link_libraries(mosflm_geometry PRIVATE hdf5_hl)
    endif()
    target_compile_definitions(mosflm_geometry PUBLIC MOSFLM_HAS_HDF5=1)
    # Vendored filter sources are C99.
    set_source_files_properties(
        ${MOSFLM_IPMOSFLM_ROOT}/lz4/lz4.c
        ${MOSFLM_IPMOSFLM_ROOT}/lz4/h5zlz4.c
        ${MOSFLM_IPMOSFLM_ROOT}/bitshuffle/bitshuffle.c
        ${MOSFLM_IPMOSFLM_ROOT}/bitshuffle/bshuf_h5filter.c
        PROPERTIES LANGUAGE C
    )
endif()

if(MOSFLM_HAS_DPS_BRAVAIS)
    set(MOSFLM_INDEXING_SOURCES
        src/indexing/BravaisClassifier.cpp
        src/indexing/BravaisSelect.cpp
        src/indexing/SpotList.cpp
        src/indexing/SpotFinder.cpp
        src/indexing/MultiFrameSpotFinder.cpp
    )
    if(MOSFLM_HAS_DPS_INDEXER)
        list(APPEND MOSFLM_INDEXING_SOURCES
            src/indexing/AutoIndexer.cpp
            src/indexing/DiffVectorIndex.cpp
            src/indexing/ControlAutoIndex.cpp
            src/indexing/ControlFindSpots.cpp
        )
    endif()
    add_library(mosflm_indexing ${MOSFLM_INDEXING_SOURCES})
    target_include_directories(mosflm_indexing
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/include
        PRIVATE
            ${MOSFLM_IPMOSFLM_ROOT}/util
            ${MOSFLM_IPMOSFLM_ROOT}/index
            ${MOSFLM_IPMOSFLM_ROOT}/include
            ${MOSFLM_IPMOSFLM_ROOT}/src/dps/index
            ${MOSFLM_IPMOSFLM_ROOT}/src/dps/util
    )
    target_link_libraries(mosflm_indexing
        PUBLIC
            mosflm_geometry
        PRIVATE
            m
    )
    # findlen.f is required by ind_gen_direct but not in libdps_index.a.
    set(MOSFLM_FINDLEN_F "${MOSFLM_IPMOSFLM_ROOT}/src/dps/index/findlen.f")
    # Order matters for single-pass ld; consumers use RESCAN (Linux) or
    # repeated archives (Apple) to resolve circular DPS symbols.
    set(MOSFLM_DPS_LINK_LIBS
        ${MOSFLM_DPS_INDEX_LIB}
        ${MOSFLM_DPS_UTIL_LIB}
    )
    if(MOSFLM_HAS_DPS_INDEXER)
        list(APPEND MOSFLM_DPS_LINK_LIBS ${MOSFLM_DPS_INDEXER_LIB})
        if(EXISTS "${MOSFLM_DPS_FFTPACK_LIB}")
            list(APPEND MOSFLM_DPS_LINK_LIBS ${MOSFLM_DPS_FFTPACK_LIB})
        endif()
        if(MOSFLM_HAS_FORTRAN AND EXISTS "${MOSFLM_FINDLEN_F}")
            target_sources(mosflm_indexing PRIVATE ${MOSFLM_FINDLEN_F})
            set_source_files_properties(${MOSFLM_FINDLEN_F}
                PROPERTIES Fortran_FORMAT FIXED
            )
            if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
                list(APPEND MOSFLM_DPS_LINK_LIBS gfortran)
            endif()
        else()
            message(WARNING
                "AutoIndexer needs findlen.f + Fortran compiler; "
                "insulin_autoindex may fail to link")
        endif()
    endif()
    # Propagate DPS archives to executables (PUBLIC) so they appear on the
    # final link line; ProcessingLab wraps them for circular resolution.
    target_link_libraries(mosflm_indexing PUBLIC ${MOSFLM_DPS_LINK_LIBS})
endif()

add_executable(mosflm_cpp_cli
    app/main.cpp
)

target_link_libraries(mosflm_cpp_cli
    PRIVATE
        mosflm_geometry
)

# Default release binary (C++ CONTROL). Distinct from F77 ``ipmosflm`` oracle.
# On Windows, CMake appends ``.exe`` → ProcessingLab.exe.
add_executable(ProcessingLab
    app/ipmosflm_main.cpp
)
target_link_libraries(ProcessingLab
    PRIVATE
        mosflm_geometry
)
if(TARGET mosflm_indexing AND MOSFLM_HAS_DPS_INDEXER)
    # Circular static deps: libdps_index ↔ libdps ↔ findlen_ in mosflm_indexing.
    # LINK_GROUP:RESCAN needs GNU ld (--start-group); Apple ld does not support
    # that feature for CXX, so repeat the archive list instead.
    if(APPLE)
        target_link_libraries(ProcessingLab PRIVATE
            mosflm_indexing
            ${MOSFLM_DPS_LINK_LIBS}
            mosflm_indexing
            ${MOSFLM_DPS_LINK_LIBS}
        )
    else()
        target_link_libraries(ProcessingLab PRIVATE
            "$<LINK_GROUP:RESCAN,mosflm_indexing,${MOSFLM_DPS_LINK_LIBS}>"
        )
    endif()
    target_compile_definitions(ProcessingLab PRIVATE MOSFLM_HAS_INDEXING=1)
endif()
set_target_properties(ProcessingLab PROPERTIES
    OUTPUT_NAME ProcessingLab
)

install(TARGETS mosflm_cpp_cli ProcessingLab
    RUNTIME DESTINATION bin
)
install(TARGETS mosflm_geometry
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
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
    # Shared assertion helpers (tests/test_util.hpp).
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tests)

    add_executable(test_unit_cell
        tests/unit/test_unit_cell.cpp
    )
    target_link_libraries(test_unit_cell PRIVATE mosflm_geometry)
    add_test(NAME unit_cell COMMAND test_unit_cell)

    add_executable(test_cell_fix
        tests/unit/test_cell_fix.cpp
    )
    target_link_libraries(test_cell_fix PRIVATE mosflm_geometry)
    add_test(NAME cell_fix COMMAND test_cell_fix)

    add_executable(test_xparm_io
        tests/unit/test_xparm_io.cpp
    )
    target_link_libraries(test_xparm_io PRIVATE mosflm_geometry)
    add_test(NAME xparm_io COMMAND test_xparm_io)

    add_executable(test_trusted_region
        tests/unit/test_trusted_region.cpp
    )
    target_link_libraries(test_trusted_region PRIVATE mosflm_geometry)
    add_test(NAME trusted_region COMMAND test_trusted_region)

    add_executable(test_detector_init
        tests/unit/test_detector_init.cpp
    )
    target_link_libraries(test_detector_init PRIVATE mosflm_geometry)
    add_test(NAME detector_init COMMAND test_detector_init)

    add_executable(test_partial_absorb
        tests/unit/test_partial_absorb.cpp
    )
    target_link_libraries(test_partial_absorb PRIVATE mosflm_geometry)
    add_test(NAME partial_absorb COMMAND test_partial_absorb)

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
    set_tests_properties(metric_invcel PROPERTIES
        ENVIRONMENT
            "MOSFLM_ANALYTICAL_FIXTURE=${CMAKE_CURRENT_SOURCE_DIR}/testdata/analytical/thaumatin_cell.txt"
    )

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

    add_executable(test_smv_image_reader
        tests/unit/test_smv_image_reader.cpp
    )
    target_link_libraries(test_smv_image_reader PRIVATE mosflm_geometry)
    add_test(
        NAME smv_image_reader
        COMMAND test_smv_image_reader
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_measurement_box
        tests/unit/test_measurement_box.cpp
    )
    target_link_libraries(test_measurement_box PRIVATE mosflm_geometry)
    add_test(NAME measurement_box COMMAND test_measurement_box)

    add_executable(test_standard_profile
        tests/unit/test_standard_profile.cpp
    )
    target_link_libraries(test_standard_profile PRIVATE mosflm_geometry)
    add_test(NAME standard_profile COMMAND test_standard_profile)

    add_executable(test_post_refinement
        tests/unit/test_post_refinement.cpp
    )
    target_link_libraries(test_post_refinement PRIVATE mosflm_geometry)
    add_test(NAME post_refinement COMMAND test_post_refinement)
    add_executable(test_prediction_divergence
        tests/unit/test_prediction_divergence.cpp
    )
    target_link_libraries(test_prediction_divergence PRIVATE mosflm_geometry)
    add_test(NAME prediction_divergence COMMAND test_prediction_divergence)

    add_executable(test_control_session
        tests/unit/test_control_session.cpp
    )
    target_link_libraries(test_control_session PRIVATE mosflm_geometry)
    add_test(
        NAME control_session
        COMMAND test_control_session
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_xinfo_xds_import
        tests/unit/test_xinfo_xds_import.cpp
    )
    target_link_libraries(test_xinfo_xds_import PRIVATE mosflm_geometry)
    add_test(
        NAME xinfo_xds_import
        COMMAND test_xinfo_xds_import
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_data_quality_stats
        tests/unit/test_data_quality_stats.cpp
    )
    target_link_libraries(test_data_quality_stats PRIVATE mosflm_geometry)
    add_test(NAME data_quality_stats COMMAND test_data_quality_stats)

    add_executable(test_image_io tests/unit/test_image_io.cpp)
    target_link_libraries(test_image_io PRIVATE mosflm_geometry)
    add_test(NAME image_io COMMAND test_image_io)

    add_executable(test_imosflm_socket tests/unit/test_imosflm_socket.cpp)
    target_link_libraries(test_imosflm_socket PRIVATE mosflm_geometry)
    add_test(NAME imosflm_socket COMMAND test_imosflm_socket)
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

    # F77 ipmosflm vs C++ ProcessPipeline on a real CBF (soft-skip if missing).
    add_executable(test_cbf_process_compare
        tests/differential/test_cbf_process_compare.cpp
    )
    target_link_libraries(test_cbf_process_compare PRIVATE mosflm_geometry)
    add_test(
        NAME cbf_process_compare
        COMMAND test_cbf_process_compare
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    set_tests_properties(cbf_process_compare PROPERTIES
        ENVIRONMENT
            "MOSFLM_IPMOSFLM=${MOSFLM_IPMOSFLM_ROOT}/bin/ipmosflm.oracle"
    )

    # F77 ipmosflm vs ProcessingLab CLI on known-indexed insulin SMV PROCESS.
    # Compares MTZ cell/SG + ASU-aware HKL recall (soft-skip if missing).
    add_executable(test_smv_process_compare
        tests/differential/test_smv_process_compare.cpp
    )
    target_link_libraries(test_smv_process_compare PRIVATE mosflm_geometry)
    add_dependencies(test_smv_process_compare ProcessingLab)
    add_test(
        NAME smv_process_compare
        COMMAND test_smv_process_compare
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    if(DEFINED ENV{CCP4} AND NOT "$ENV{CCP4}" STREQUAL "")
        set(_MOSFLM_SMV_CCP4 "$ENV{CCP4}")
    else()
        set(_MOSFLM_SMV_CCP4 "/Applications/ccp4-9")
    endif()
    set_tests_properties(smv_process_compare PROPERTIES
        ENVIRONMENT
            "MOSFLM_PROCESSINGLAB=${CMAKE_CURRENT_BINARY_DIR}/ProcessingLab${CMAKE_EXECUTABLE_SUFFIX};MOSFLM_IPMOSFLM=${MOSFLM_IPMOSFLM_ROOT}/bin/ipmosflm.oracle;CCP4=${_MOSFLM_SMV_CCP4}"
        TIMEOUT 300
    )

    # Thaumatin PL process MTZ vs DIALS integrated.mtz (soft-skip if missing).
    add_executable(test_thaumatin_dials_mtz_compare
        tests/differential/test_thaumatin_dials_mtz_compare.cpp
    )
    target_link_libraries(test_thaumatin_dials_mtz_compare PRIVATE mosflm_geometry)
    add_test(
        NAME thaumatin_dials_mtz_compare
        COMMAND test_thaumatin_dials_mtz_compare
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
    )
    set_tests_properties(thaumatin_dials_mtz_compare PROPERTIES TIMEOUT 120)

    # ProcessingLab vs DIALS on known-indexed insulin SMV (soft-skip if no DIALS).
    add_executable(test_smv_dials_compare
        tests/differential/test_smv_dials_compare.cpp
    )
    target_link_libraries(test_smv_dials_compare PRIVATE mosflm_geometry)
    add_dependencies(test_smv_dials_compare ProcessingLab)
    add_test(
        NAME smv_dials_compare
        COMMAND test_smv_dials_compare
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    set_tests_properties(smv_dials_compare PROPERTIES
        ENVIRONMENT
            "MOSFLM_PROCESSINGLAB=${CMAKE_CURRENT_BINARY_DIR}/ProcessingLab${CMAKE_EXECUTABLE_SUFFIX};MOSFLM_DIALS_BIN=${_MOSFLM_SMV_CCP4}/bin;CCP4=${_MOSFLM_SMV_CCP4}"
        TIMEOUT 600
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

    add_executable(test_lab_frame
        tests/unit/test_lab_frame.cpp
    )
    target_link_libraries(test_lab_frame PRIVATE mosflm_geometry)
    add_test(NAME lab_frame COMMAND test_lab_frame)

    if(MOSFLM_FORTRAN_TESTS)
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

        add_executable(test_predict_gendata_pack
            fortran/predict_cpp_bridge.f90
            fortran/predict_gendata_pack.f90
            fortran/test_predict_gendata_pack.f90
        )
        target_link_libraries(test_predict_gendata_pack PRIVATE mosflm_geometry)
        set_target_properties(test_predict_gendata_pack PROPERTIES
            LINKER_LANGUAGE CXX
            Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_gendata_mods
        )
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
            target_link_libraries(test_predict_gendata_pack PRIVATE gfortran)
        endif()
        add_test(NAME predict_gendata_pack COMMAND test_predict_gendata_pack)

        add_executable(test_predict_cpp_hook
            fortran/predict_cpp_bridge.f90
            fortran/predict_gendata_pack.f90
            fortran/predict_cpp_hook.f90
            fortran/test_predict_cpp_hook.f90
        )
        target_link_libraries(test_predict_cpp_hook PRIVATE mosflm_geometry)
        set_target_properties(test_predict_cpp_hook PROPERTIES
            LINKER_LANGUAGE CXX
            Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_hook_mods
        )
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
            target_link_libraries(test_predict_cpp_hook PRIVATE gfortran)
        endif()
        add_test(NAME predict_cpp_hook COMMAND test_predict_cpp_hook)

        add_executable(test_predict_cpp_hook_f77
            fortran/predict_cpp_bridge.f90
            fortran/predict_gendata_pack.f90
            fortran/predict_cpp_hook.f90
            fortran/test_predict_cpp_hook_f77.f
        )
        set_source_files_properties(fortran/test_predict_cpp_hook_f77.f
            PROPERTIES Fortran_FORMAT FIXED
        )
        target_link_libraries(test_predict_cpp_hook_f77 PRIVATE mosflm_geometry)
        set_target_properties(test_predict_cpp_hook_f77 PROPERTIES
            LINKER_LANGUAGE CXX
            Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_hook_f77_mods
        )
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
            target_link_libraries(test_predict_cpp_hook_f77 PRIVATE gfortran)
        endif()
        add_test(NAME predict_cpp_hook_f77 COMMAND test_predict_cpp_hook_f77)

        add_executable(test_predict_cpp_from_common
            fortran/predict_cpp_bridge.f90
            fortran/predict_gendata_pack.f90
            fortran/predict_cpp_hook.f90
            fortran/predict_cpp_from_common.f90
            fortran/test_predict_cpp_from_common.f90
        )
        target_link_libraries(test_predict_cpp_from_common PRIVATE mosflm_geometry)
        set_target_properties(test_predict_cpp_from_common PROPERTIES
            LINKER_LANGUAGE CXX
            Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_from_common_mods
        )
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
            target_link_libraries(test_predict_cpp_from_common PRIVATE gfortran)
        endif()
        add_test(NAME predict_cpp_from_common COMMAND test_predict_cpp_from_common)
    endif()

    # Amalgamation splice library — opt-in (Phase D / D6). Default builds do
    # not produce artifacts that link mosflm_all_ip_inc_noX.o.
    if(MOSFLM_HAS_FORTRAN AND MOSFLM_FORTRAN_ORACLE)
        # Static archive for Autotools
        # ``make MOSFLM_FORTRAN_ORACLE=1 MOSFLM_CPP_PREDICT=1 ipmosflm``.
        add_library(mosflm_predict_hook STATIC
            fortran/predict_cpp_bridge.f90
            fortran/predict_gendata_pack.f90
            fortran/predict_cpp_hook.f90
            fortran/predict_cpp_from_common.f90
        )
        target_link_libraries(mosflm_predict_hook PUBLIC mosflm_geometry)
        set_target_properties(mosflm_predict_hook PROPERTIES
            Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_hook_lib_mods
            POSITION_INDEPENDENT_CODE ON
        )
        message(STATUS "MOSFLM_FORTRAN_ORACLE=ON: building mosflm_predict_hook")
    elseif(MOSFLM_FORTRAN_ORACLE AND NOT MOSFLM_HAS_FORTRAN)
        message(WARNING
            "MOSFLM_FORTRAN_ORACLE=ON but no Fortran compiler; "
            "skipping mosflm_predict_hook")
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

    add_executable(test_insulin_hkl_overlap
        tests/differential/test_insulin_hkl_overlap.cpp
    )
    target_link_libraries(test_insulin_hkl_overlap PRIVATE mosflm_geometry)
    add_test(
        NAME insulin_hkl_overlap
        COMMAND test_insulin_hkl_overlap
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_insulin_xdet
        tests/differential/test_insulin_xdet.cpp
    )
    target_link_libraries(test_insulin_xdet PRIVATE mosflm_geometry)
    add_test(
        NAME insulin_xdet
        COMMAND test_insulin_xdet
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_insulin_process_mtz
        tests/differential/test_insulin_process_mtz.cpp
    )
    target_link_libraries(test_insulin_process_mtz PRIVATE mosflm_geometry)
    add_test(
        NAME insulin_process_mtz
        COMMAND test_insulin_process_mtz
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_insulin_cpp_process
        tests/differential/test_insulin_cpp_process.cpp
    )
    target_link_libraries(test_insulin_cpp_process PRIVATE mosflm_geometry)
    add_test(
        NAME insulin_cpp_process
        COMMAND test_insulin_cpp_process
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_insulin_isi
        tests/differential/test_insulin_isi.cpp
    )
    target_link_libraries(test_insulin_isi PRIVATE mosflm_geometry)
    add_test(
        NAME insulin_isi
        COMMAND test_insulin_isi
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_insulin_profile
        tests/differential/test_insulin_profile.cpp
    )
    target_link_libraries(test_insulin_profile PRIVATE mosflm_geometry)
    add_test(
        NAME insulin_profile
        COMMAND test_insulin_profile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_insulin_postref
        tests/differential/test_insulin_postref.cpp
    )
    target_link_libraries(test_insulin_postref PRIVATE mosflm_geometry)
    add_test(
        NAME insulin_postref
        COMMAND test_insulin_postref
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_insulin_control
        tests/differential/test_insulin_control.cpp
    )
    target_link_libraries(test_insulin_control PRIVATE mosflm_geometry)
    add_test(
        NAME insulin_control
        COMMAND test_insulin_control
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_integrationlab
        tests/differential/test_ipmosflm_thin.cpp
    )
    add_dependencies(test_integrationlab ProcessingLab)
    add_test(
        NAME integrationlab
        COMMAND test_integrationlab
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    set_tests_properties(integrationlab PROPERTIES
        ENVIRONMENT
            "MOSFLM_PROCESSINGLAB=${CMAKE_CURRENT_BINARY_DIR}/ProcessingLab${CMAKE_EXECUTABLE_SUFFIX}"
    )

    # Wall-time: ProcessingLab (C++) vs F77 ipmosflm on known-indexed PROCESS.
    # Soft-skips inside the binary if either side is missing / wrong arch.
    add_executable(test_ipmosflm_perf
        tests/differential/test_ipmosflm_perf.cpp
    )
    add_dependencies(test_ipmosflm_perf ProcessingLab)
    add_test(
        NAME ipmosflm_perf
        COMMAND test_ipmosflm_perf
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    if(DEFINED ENV{CCP4} AND NOT "$ENV{CCP4}" STREQUAL "")
        set(_MOSFLM_PERF_CCP4 "$ENV{CCP4}")
    else()
        set(_MOSFLM_PERF_CCP4 "/apps/ccp4/ccp4-9")
    endif()
    set_tests_properties(ipmosflm_perf PROPERTIES
        ENVIRONMENT
            "MOSFLM_PROCESSINGLAB=${CMAKE_CURRENT_BINARY_DIR}/ProcessingLab${CMAKE_EXECUTABLE_SUFFIX};MOSFLM_IPMOSFLM=${_MOSFLM_PERF_CCP4}/bin/ipmosflm;CCP4=${_MOSFLM_PERF_CCP4}"
        TIMEOUT 600
    )

    if(MOSFLM_HAS_DPS_BRAVAIS)
        add_executable(test_bravais_classifier
            tests/differential/test_bravais_classifier.cpp
        )
        target_link_libraries(test_bravais_classifier PRIVATE mosflm_indexing)
        add_test(NAME bravais_classifier COMMAND test_bravais_classifier)

        add_executable(test_bravais_select
            tests/unit/test_bravais_select.cpp
        )
        target_link_libraries(test_bravais_select PRIVATE mosflm_indexing)
        add_test(
            NAME bravais_select
            COMMAND test_bravais_select
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_executable(test_spot_list
            tests/unit/test_spot_list.cpp
        )
        target_link_libraries(test_spot_list PRIVATE mosflm_indexing)
        add_test(
            NAME spot_list
            COMMAND test_spot_list
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_executable(test_spot_finder
            tests/unit/test_spot_finder.cpp
        )
        target_link_libraries(test_spot_finder PRIVATE mosflm_indexing)
        add_test(
            NAME spot_finder
            COMMAND test_spot_finder
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
    endif()

    if(MOSFLM_HAS_DPS_INDEXER)
        add_executable(test_insulin_autoindex
            tests/differential/test_insulin_autoindex.cpp
        )
        target_link_libraries(test_insulin_autoindex PRIVATE mosflm_indexing)
        add_test(
            NAME insulin_autoindex
            COMMAND test_insulin_autoindex
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_executable(test_control_autoindex
            tests/differential/test_control_autoindex.cpp
        )
        target_link_libraries(test_control_autoindex PRIVATE mosflm_indexing)
        add_test(
            NAME control_autoindex
            COMMAND test_control_autoindex
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_executable(test_diff_vector_index
            tests/unit/test_diff_vector_index.cpp
        )
        target_link_libraries(test_diff_vector_index PRIVATE mosflm_indexing)
        add_test(
            NAME diff_vector_index
            COMMAND test_diff_vector_index
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_executable(test_76c1b_autoindex
            tests/differential/test_76c1b_autoindex.cpp
        )
        target_link_libraries(test_76c1b_autoindex PRIVATE mosflm_indexing)
        add_test(
            NAME 76c1b_autoindex
            COMMAND test_76c1b_autoindex
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
    endif()
endif()
