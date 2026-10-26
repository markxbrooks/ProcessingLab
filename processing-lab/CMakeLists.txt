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
    "Build processinglab_predict_hook for Autotools amalgamation splice (opt-in; D6)"
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

# Prefer a system / Homebrew Eigen (header-only). GitLab git clone is a
# last resort and often fails under load; use a tarball URL instead of git.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(Eigen3_FOUND)
    message(STATUS "Using system Eigen ${Eigen3_VERSION}")
else()
    FetchContent_Declare(
        eigen
        URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
        URL_HASH SHA256=8586084f71f9bde545ee7fa6d97872b0c51dd0edc0005e9ba67b37c35a58089e
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING_SAVED "${BUILD_TESTING}")
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(eigen)
    set(BUILD_TESTING "${BUILD_TESTING_SAVED}" CACHE BOOL "" FORCE)
endif()

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

# Vendored CBFlib (Option A) — sibling of processing-lab under the repo root.
# Stale IDE caches from the mosflm-cpp → processing-lab rename leave
# MOSFLM_IPMOSFLM_ROOT as ".../mosflm-cpp/.."; CMake EXISTS rejects that
# unresolved path even though libcbf.a is present after normalization.
get_filename_component(_mosflm_repo_root "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
if(NOT DEFINED CACHE{MOSFLM_IPMOSFLM_ROOT}
   OR "${MOSFLM_IPMOSFLM_ROOT}" MATCHES "mosflm-cpp"
   OR NOT EXISTS "${MOSFLM_IPMOSFLM_ROOT}/cbf/lib/libcbf.a")
    set(MOSFLM_IPMOSFLM_ROOT "${_mosflm_repo_root}" CACHE PATH
        "ProcessingLab repository root (CBF + DPS libs)" FORCE)
endif()
get_filename_component(MOSFLM_IPMOSFLM_ROOT "${MOSFLM_IPMOSFLM_ROOT}" ABSOLUTE)
set(MOSFLM_CBF_INCLUDE_DIR "${MOSFLM_IPMOSFLM_ROOT}/cbf/include")
set(MOSFLM_CBF_LIB "${MOSFLM_IPMOSFLM_ROOT}/cbf/lib/libcbf.a")
if(NOT EXISTS "${MOSFLM_CBF_LIB}")
    message(FATAL_ERROR
        "CBFlib not found at ${MOSFLM_CBF_LIB}\n"
        "Build it first: make -C ${MOSFLM_IPMOSFLM_ROOT}/cbf  (or: make cbf-lib)")
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

# Optional OpenMP for frame-level FINDSPOTS / PROCESS parallelism.
# Apple Clang needs Homebrew libomp (-Xpreprocessor -fopenmp); plain
# find_package(OpenMP) often fails without hints.
if(APPLE AND NOT OpenMP_CXX_FOUND)
    foreach(_omp_prefix
            "$ENV{HOMEBREW_PREFIX}/opt/libomp"
            /opt/homebrew/opt/libomp
            /usr/local/opt/libomp)
        if(EXISTS "${_omp_prefix}/lib/libomp.dylib")
            set(OpenMP_CXX_FLAGS
                "-Xpreprocessor -fopenmp -I${_omp_prefix}/include"
                CACHE STRING "OpenMP CXX flags" FORCE)
            set(OpenMP_C_FLAGS
                "-Xpreprocessor -fopenmp -I${_omp_prefix}/include"
                CACHE STRING "OpenMP C flags" FORCE)
            set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "OpenMP CXX libs" FORCE)
            set(OpenMP_C_LIB_NAMES "omp" CACHE STRING "OpenMP C libs" FORCE)
            set(OpenMP_omp_LIBRARY "${_omp_prefix}/lib/libomp.dylib"
                CACHE FILEPATH "libomp" FORCE)
            break()
        endif()
    endforeach()
endif()
find_package(OpenMP)
set(MOSFLM_HAS_OPENMP OFF)
if(OpenMP_CXX_FOUND)
    set(MOSFLM_HAS_OPENMP ON)
    message(STATUS "OpenMP found — enabling NPROC frame parallelism")
else()
    message(WARNING
        "OpenMP not found — NPROC runs serial "
        "(macOS: brew install libomp)")
endif()

# Optional FFTW for the indexing Patterson FFT.
# Preferred path: Autotools ./configure detects FFTW and injects flags via
# Makefile (-DMOSFLM_FFTW_FROM_CONFIGURE=ON). Direct CMake find_* is legacy.
option(MOSFLM_FFTW_FROM_CONFIGURE
    "FFTW on/off already decided by Autotools configure"
    OFF
)
option(MOSFLM_HAS_FFTW "Enable FFTW3 for indexing FFTs" OFF)
set(MOSFLM_AUTOTOOLS_FFTW OFF)
if(MOSFLM_FFTW_FROM_CONFIGURE)
    if(MOSFLM_HAS_FFTW)
        set(MOSFLM_AUTOTOOLS_FFTW ON)
        message(STATUS
            "FFTW from Autotools — indexing uses fftw_plan_dft_3d "
            "(CFLAGS='${MOSFLM_FFTW_CFLAGS}' LIBS='${MOSFLM_FFTW_LIBS}')")
    else()
        message(STATUS
            "FFTW disabled by Autotools configure — built-in radix-2 FFT")
    endif()
else()
    find_package(PkgConfig QUIET)
    set(MOSFLM_HAS_FFTW OFF)
    if(PkgConfig_FOUND)
        pkg_check_modules(FFTW3 QUIET IMPORTED_TARGET fftw3)
    endif()
    if(NOT TARGET PkgConfig::FFTW3)
        find_path(FFTW3_INCLUDE_DIR fftw3.h
            HINTS
                "$ENV{HOMEBREW_PREFIX}/include"
                /opt/homebrew/include
                /usr/local/include
        )
        find_library(FFTW3_LIBRARY NAMES fftw3
            HINTS
                "$ENV{HOMEBREW_PREFIX}/lib"
                /opt/homebrew/lib
                /usr/local/lib
        )
    endif()
    if(TARGET PkgConfig::FFTW3 OR (FFTW3_INCLUDE_DIR AND FFTW3_LIBRARY))
        set(MOSFLM_HAS_FFTW ON)
        message(STATUS
            "FFTW found (legacy CMake search) — indexing uses fftw_plan_dft_3d")
    else()
        message(STATUS
            "FFTW not found — indexing uses the built-in radix-2 FFT "
            "(./configure --with-fftw, or: apt install libfftw3-dev / brew install fftw)")
    endif()
endif()

set(MOSFLM_GEOMETRY_SOURCES
    src/geometry/UnitCell.cpp
    src/geometry/CellFix.cpp
    src/geometry/MetricTensorProjector.cpp
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
    src/diffraction/PartialSummation.cpp
    src/diffraction/MeasurementBox.cpp
    src/diffraction/TrustedRegion.cpp
    src/diffraction/Defpix.cpp
    src/diffraction/Strategy.cpp
    src/diffraction/DetectorInit.cpp
    src/diffraction/StandardProfile.cpp
    src/diffraction/Shoebox.cpp
    src/diffraction/Profile3d.cpp
    src/diffraction/ProfileFitter.cpp
    src/diffraction/LpCorrection.cpp
    src/diffraction/DataQualityStats.cpp
    src/diffraction/SystematicAbsences.cpp
    src/refinement/PostRefinement.cpp
    src/refinement/PositionalRefinement.cpp
    src/control/ControlSession.cpp
    src/io/ImportedSweepGeometry.cpp
    src/io/XinfoImporter.cpp
    src/io/ExperimentModel.cpp
    src/io/Json.cpp
    src/io/RepeatInputs.cpp
    src/io/RunRecord.cpp
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

add_library(processinglab_geometry ${MOSFLM_GEOMETRY_SOURCES})

target_include_directories(processinglab_geometry
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${gemmi_SOURCE_DIR}/include
    PRIVATE
        ${MOSFLM_CBF_INCLUDE_DIR}
)

target_link_libraries(processinglab_geometry
    PUBLIC
        Eigen3::Eigen
        ZLIB::ZLIB
    PRIVATE
        ${MOSFLM_CBF_LIB}
        m
)

target_compile_definitions(processinglab_geometry
    PUBLIC
        GEMMI_WRITE_IMPLEMENTATION
)

if(MOSFLM_HAS_OPENMP)
    target_link_libraries(processinglab_geometry PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(processinglab_geometry
        PUBLIC PROCESSINGLAB_HAS_OPENMP=1
    )
endif()

if(MOSFLM_HAS_HDF5)
    target_include_directories(processinglab_geometry PRIVATE
        ${HDF5_INCLUDE_DIRS}
        ${MOSFLM_IPMOSFLM_ROOT}/lz4
        ${MOSFLM_IPMOSFLM_ROOT}/bitshuffle
    )
    # ${HDF5_LIBRARIES} is often only libhdf5; H5LT lives in libhdf5_hl.
    target_link_libraries(processinglab_geometry PRIVATE ${HDF5_LIBRARIES})
    if(TARGET hdf5::hdf5_hl)
        target_link_libraries(processinglab_geometry PRIVATE hdf5::hdf5_hl)
    elseif(HDF5_HL_LIBRARIES)
        target_link_libraries(processinglab_geometry PRIVATE ${HDF5_HL_LIBRARIES})
    else()
        target_link_libraries(processinglab_geometry PRIVATE hdf5_hl)
    endif()
    target_compile_definitions(processinglab_geometry PUBLIC MOSFLM_HAS_HDF5=1)
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
            src/indexing/Fft3dIndex.cpp
            src/indexing/LocalIndexAssign.cpp
            src/indexing/MaxCell.cpp
            src/indexing/IdxrefRefine.cpp
            src/indexing/SpotToRlp.cpp
            src/indexing/SupercellReduce.cpp
            src/indexing/ControlAutoIndex.cpp
            src/indexing/ControlFindSpots.cpp
        )
    endif()
    add_library(processinglab_indexing ${MOSFLM_INDEXING_SOURCES})
    target_include_directories(processinglab_indexing
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/include
        PRIVATE
            ${MOSFLM_IPMOSFLM_ROOT}/util
            ${MOSFLM_IPMOSFLM_ROOT}/index
            ${MOSFLM_IPMOSFLM_ROOT}/include
            ${MOSFLM_IPMOSFLM_ROOT}/src/dps/index
            ${MOSFLM_IPMOSFLM_ROOT}/src/dps/util
    )
    target_link_libraries(processinglab_indexing
        PUBLIC
            processinglab_geometry
        PRIVATE
            m
    )
    if(MOSFLM_HAS_FFTW)
        if(MOSFLM_AUTOTOOLS_FFTW)
            if(MOSFLM_FFTW_CFLAGS)
                separate_arguments(_mosflm_fftw_cflags NATIVE_COMMAND
                    "${MOSFLM_FFTW_CFLAGS}")
                target_compile_options(processinglab_indexing PRIVATE
                    ${_mosflm_fftw_cflags})
            endif()
            if(MOSFLM_FFTW_LIBS)
                separate_arguments(_mosflm_fftw_libs NATIVE_COMMAND
                    "${MOSFLM_FFTW_LIBS}")
                target_link_libraries(processinglab_indexing PRIVATE
                    ${_mosflm_fftw_libs})
            else()
                target_link_libraries(processinglab_indexing PRIVATE fftw3)
            endif()
        elseif(TARGET PkgConfig::FFTW3)
            target_link_libraries(processinglab_indexing PRIVATE PkgConfig::FFTW3)
        else()
            target_include_directories(processinglab_indexing PRIVATE
                ${FFTW3_INCLUDE_DIR})
            target_link_libraries(processinglab_indexing PRIVATE ${FFTW3_LIBRARY})
        endif()
        target_compile_definitions(processinglab_indexing PRIVATE MOSFLM_HAS_FFTW=1)
    endif()
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
            target_sources(processinglab_indexing PRIVATE ${MOSFLM_FINDLEN_F})
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
    target_link_libraries(processinglab_indexing PUBLIC ${MOSFLM_DPS_LINK_LIBS})
endif()

add_executable(processinglab_dev
    app/main.cpp
)

target_link_libraries(processinglab_dev
    PRIVATE
        processinglab_geometry
)

# Default release binary (C++ CONTROL). Distinct from F77 ``ipmosflm`` oracle.
# On Windows, CMake appends ``.exe`` → ProcessingLab.exe.
add_executable(ProcessingLab
    app/ipmosflm_main.cpp
)
target_link_libraries(ProcessingLab
    PRIVATE
        processinglab_geometry
)
if(TARGET processinglab_indexing AND MOSFLM_HAS_DPS_INDEXER)
    # Circular static deps: libdps_index ↔ libdps ↔ findlen_ in processinglab_indexing.
    # LINK_GROUP:RESCAN needs GNU ld (--start-group); Apple ld does not support
    # that feature for CXX, so repeat the archive list instead.
    if(APPLE)
        target_link_libraries(ProcessingLab PRIVATE
            processinglab_indexing
            ${MOSFLM_DPS_LINK_LIBS}
            processinglab_indexing
            ${MOSFLM_DPS_LINK_LIBS}
        )
    else()
        target_link_libraries(ProcessingLab PRIVATE
            "$<LINK_GROUP:RESCAN,processinglab_indexing,${MOSFLM_DPS_LINK_LIBS}>"
        )
    endif()
    target_compile_definitions(ProcessingLab PRIVATE MOSFLM_HAS_INDEXING=1)
endif()
set_target_properties(ProcessingLab PROPERTIES
    OUTPUT_NAME ProcessingLab
)

install(TARGETS ProcessingLab
    RUNTIME DESTINATION bin
)
install(TARGETS processinglab_geometry
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
    # Test outputs go under the build tree (mosflm_test::output_dir), never
    # into source, testdata or image directories.
    add_compile_definitions(PL_TEST_OUTPUT_ROOT="${CMAKE_CURRENT_BINARY_DIR}/test_output")

    add_executable(test_unit_cell
        tests/unit/test_unit_cell.cpp
    )
    target_link_libraries(test_unit_cell PRIVATE processinglab_geometry)
    add_test(NAME unit_cell COMMAND test_unit_cell)

    add_executable(test_cell_fix
        tests/unit/test_cell_fix.cpp
    )
    target_link_libraries(test_cell_fix PRIVATE processinglab_geometry)
    add_test(NAME cell_fix COMMAND test_cell_fix)

    add_executable(test_metric_tensor_projector
        tests/unit/test_metric_tensor_projector.cpp
    )
    target_link_libraries(test_metric_tensor_projector PRIVATE processinglab_geometry)
    add_test(NAME metric_tensor_projector COMMAND test_metric_tensor_projector)

    add_executable(test_rotation_lp
        tests/unit/test_rotation_lp.cpp
    )
    target_link_libraries(test_rotation_lp PRIVATE processinglab_geometry)
    add_test(NAME rotation_lp COMMAND test_rotation_lp)

    add_executable(test_positional_refinement
        tests/unit/test_positional_refinement.cpp
    )
    target_link_libraries(test_positional_refinement PRIVATE processinglab_geometry)
    add_test(NAME positional_refinement COMMAND test_positional_refinement)

    add_executable(test_partial_summation
        tests/unit/test_partial_summation.cpp
    )
    target_link_libraries(test_partial_summation PRIVATE processinglab_geometry)
    add_test(NAME partial_summation COMMAND test_partial_summation)

    add_executable(test_xparm_io
        tests/unit/test_xparm_io.cpp
    )
    target_link_libraries(test_xparm_io PRIVATE processinglab_geometry)
    add_test(NAME xparm_io COMMAND test_xparm_io)

    add_executable(test_trusted_region
        tests/unit/test_trusted_region.cpp
    )
    target_link_libraries(test_trusted_region PRIVATE processinglab_geometry)
    add_test(NAME trusted_region COMMAND test_trusted_region)

    add_executable(test_defpix
        tests/unit/test_defpix.cpp
    )
    target_link_libraries(test_defpix PRIVATE processinglab_geometry)
    add_test(NAME defpix COMMAND test_defpix)

    add_executable(test_strategy
        tests/unit/test_strategy.cpp
    )
    target_link_libraries(test_strategy PRIVATE processinglab_geometry)
    add_test(NAME strategy COMMAND test_strategy)

    add_executable(test_detector_init
        tests/unit/test_detector_init.cpp
    )
    target_link_libraries(test_detector_init PRIVATE processinglab_geometry)
    add_test(NAME detector_init COMMAND test_detector_init)

    add_executable(test_partial_absorb
        tests/unit/test_partial_absorb.cpp
    )
    target_link_libraries(test_partial_absorb PRIVATE processinglab_geometry)
    add_test(NAME partial_absorb COMMAND test_partial_absorb)

    add_executable(test_orientation
        tests/unit/test_orientation.cpp
    )
    target_link_libraries(test_orientation PRIVATE processinglab_geometry)
    add_test(NAME orientation COMMAND test_orientation)

    add_executable(test_detector
        tests/unit/test_detector.cpp
    )
    target_link_libraries(test_detector PRIVATE processinglab_geometry)
    add_test(NAME detector COMMAND test_detector)

    add_executable(test_metric_invcel
        tests/differential/test_metric_invcel.cpp
    )
    target_link_libraries(test_metric_invcel PRIVATE processinglab_geometry)
    add_test(NAME metric_invcel COMMAND test_metric_invcel)
    set_tests_properties(metric_invcel PROPERTIES
        ENVIRONMENT
            "MOSFLM_ANALYTICAL_FIXTURE=${CMAKE_CURRENT_SOURCE_DIR}/testdata/analytical/thaumatin_cell.txt"
    )

    add_executable(test_gemmi_crystal
        tests/differential/test_gemmi_crystal.cpp
    )
    target_link_libraries(test_gemmi_crystal PRIVATE processinglab_geometry)
    add_test(NAME gemmi_crystal COMMAND test_gemmi_crystal)

    add_executable(test_cbf_image_reader
        tests/unit/test_cbf_image_reader.cpp
    )
    target_link_libraries(test_cbf_image_reader PRIVATE processinglab_geometry)
    add_test(
        NAME cbf_image_reader
        COMMAND test_cbf_image_reader
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_smv_image_reader
        tests/unit/test_smv_image_reader.cpp
    )
    target_link_libraries(test_smv_image_reader PRIVATE processinglab_geometry)
    add_test(
        NAME smv_image_reader
        COMMAND test_smv_image_reader
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_measurement_box
        tests/unit/test_measurement_box.cpp
    )
    target_link_libraries(test_measurement_box PRIVATE processinglab_geometry)
    add_test(NAME measurement_box COMMAND test_measurement_box)

    add_executable(test_standard_profile
        tests/unit/test_standard_profile.cpp
    )
    target_link_libraries(test_standard_profile PRIVATE processinglab_geometry)
    add_test(NAME standard_profile COMMAND test_standard_profile)

    add_executable(test_profile3d
        tests/unit/test_profile3d.cpp
    )
    target_link_libraries(test_profile3d PRIVATE processinglab_geometry)
    add_test(NAME profile3d COMMAND test_profile3d)

    add_executable(test_post_refinement
        tests/unit/test_post_refinement.cpp
    )
    target_link_libraries(test_post_refinement PRIVATE processinglab_geometry)
    add_test(NAME post_refinement COMMAND test_post_refinement)
    add_executable(test_prediction_divergence
        tests/unit/test_prediction_divergence.cpp
    )
    target_link_libraries(test_prediction_divergence PRIVATE processinglab_geometry)
    add_test(NAME prediction_divergence COMMAND test_prediction_divergence)

    add_executable(test_control_session
        tests/unit/test_control_session.cpp
    )
    target_link_libraries(test_control_session PRIVATE processinglab_geometry)
    add_test(
        NAME control_session
        COMMAND test_control_session
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_openmp_nproc
        tests/unit/test_openmp_nproc.cpp
    )
    target_link_libraries(test_openmp_nproc PRIVATE processinglab_geometry)
    add_test(NAME openmp_nproc COMMAND test_openmp_nproc)

    add_executable(test_xinfo_xds_import
        tests/unit/test_xinfo_xds_import.cpp
    )
    target_link_libraries(test_xinfo_xds_import PRIVATE processinglab_geometry)
    add_test(
        NAME xinfo_xds_import
        COMMAND test_xinfo_xds_import
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_run_record
        tests/unit/test_run_record.cpp
    )
    target_link_libraries(test_run_record PRIVATE processinglab_geometry)
    add_test(
        NAME run_record
        COMMAND test_run_record
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_experiment_model
        tests/unit/test_experiment_model.cpp
    )
    target_link_libraries(test_experiment_model PRIVATE processinglab_geometry)
    add_test(
        NAME experiment_model
        COMMAND test_experiment_model
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_data_quality_stats
        tests/unit/test_data_quality_stats.cpp
    )
    target_link_libraries(test_data_quality_stats PRIVATE processinglab_geometry)
    add_test(NAME data_quality_stats COMMAND test_data_quality_stats)

    add_executable(test_image_io tests/unit/test_image_io.cpp)
    target_link_libraries(test_image_io PRIVATE processinglab_geometry)
    add_test(NAME image_io COMMAND test_image_io)

    if(MOSFLM_HAS_HDF5)
        # Writes synthetic NeXus files, so it needs HDF5 and bitshuffle directly.
        add_executable(test_hdf5_reader tests/unit/test_hdf5_reader.cpp)
        target_include_directories(test_hdf5_reader PRIVATE
            ${HDF5_INCLUDE_DIRS}
            ${MOSFLM_IPMOSFLM_ROOT}/bitshuffle
        )
        target_link_libraries(test_hdf5_reader PRIVATE processinglab_geometry
            ${HDF5_LIBRARIES})
        if(TARGET hdf5::hdf5_hl)
            target_link_libraries(test_hdf5_reader PRIVATE hdf5::hdf5_hl)
        elseif(HDF5_HL_LIBRARIES)
            target_link_libraries(test_hdf5_reader PRIVATE ${HDF5_HL_LIBRARIES})
        else()
            target_link_libraries(test_hdf5_reader PRIVATE hdf5_hl)
        endif()
        add_test(NAME hdf5_reader COMMAND test_hdf5_reader)
    endif()

    add_executable(test_imosflm_socket tests/unit/test_imosflm_socket.cpp)
    target_link_libraries(test_imosflm_socket PRIVATE processinglab_geometry)
    add_test(NAME imosflm_socket COMMAND test_imosflm_socket)
    set_tests_properties(cbf_image_reader PROPERTIES
        ENVIRONMENT
            "MOSFLM_CBF_FIXTURE=${CMAKE_CURRENT_SOURCE_DIR}/testdata/cbf/sample_minicbf.cbf"
    )

    add_executable(test_predict_stub
        tests/differential/test_predict_stub.cpp
    )
    target_link_libraries(test_predict_stub PRIVATE processinglab_geometry)
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
    target_link_libraries(test_summation_integrate PRIVATE processinglab_geometry)
    add_test(
        NAME summation_integrate
        COMMAND test_summation_integrate
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_dstar_rotmat
        tests/unit/test_dstar_rotmat.cpp
    )
    target_link_libraries(test_dstar_rotmat PRIVATE processinglab_geometry)
    add_test(NAME dstar_rotmat COMMAND test_dstar_rotmat)

    add_executable(test_clcalc_matset
        tests/unit/test_clcalc_matset.cpp
    )
    target_link_libraries(test_clcalc_matset PRIVATE processinglab_geometry)
    add_test(NAME clcalc_matset COMMAND test_clcalc_matset)

    add_executable(test_vertical_slice
        tests/differential/test_vertical_slice.cpp
    )
    target_link_libraries(test_vertical_slice PRIVATE processinglab_geometry)
    add_test(NAME vertical_slice COMMAND test_vertical_slice)

    add_executable(test_predict_osc
        tests/differential/test_predict_osc.cpp
    )
    target_link_libraries(test_predict_osc PRIVATE processinglab_geometry)
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
    target_link_libraries(test_cbf_vertical_slice PRIVATE processinglab_geometry)
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
    target_link_libraries(test_cbf_process_compare PRIVATE processinglab_geometry)
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
    # ProcessingLab input modes (run / imports / script) on CCP4 insulin.
    add_executable(test_run_modes
        tests/differential/test_run_modes.cpp
    )
    target_link_libraries(test_run_modes PRIVATE processinglab_geometry)
    target_compile_definitions(test_run_modes PRIVATE
        PROCESSINGLAB_BIN="$<TARGET_FILE:ProcessingLab>"
    )
    add_dependencies(test_run_modes ProcessingLab)
    add_test(
        NAME run_modes
        COMMAND test_run_modes
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_smv_process_compare
        tests/differential/test_smv_process_compare.cpp
    )
    target_link_libraries(test_smv_process_compare PRIVATE processinglab_geometry)
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
    target_link_libraries(test_thaumatin_dials_mtz_compare PRIVATE processinglab_geometry)
    add_test(
        NAME thaumatin_dials_mtz_compare
        COMMAND test_thaumatin_dials_mtz_compare
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
    )
    set_tests_properties(thaumatin_dials_mtz_compare PROPERTIES TIMEOUT 120)

    add_executable(test_thaumatin_profile3d
        tests/differential/test_thaumatin_profile3d.cpp
    )
    target_link_libraries(test_thaumatin_profile3d PRIVATE processinglab_geometry)
    add_test(
        NAME thaumatin_profile3d
        COMMAND test_thaumatin_profile3d
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(thaumatin_profile3d PROPERTIES TIMEOUT 600)

    # ProcessingLab vs DIALS on known-indexed insulin SMV (soft-skip if no DIALS).
    add_executable(test_smv_dials_compare
        tests/differential/test_smv_dials_compare.cpp
    )
    target_link_libraries(test_smv_dials_compare PRIVATE processinglab_geometry)
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
    target_link_libraries(test_profile_fitter PRIVATE processinglab_geometry)
    add_test(NAME profile_fitter COMMAND test_profile_fitter)

    add_executable(test_c_api_lp
        tests/unit/test_c_api_lp.cpp
    )
    target_link_libraries(test_c_api_lp PRIVATE processinglab_geometry)
    add_test(NAME c_api_lp COMMAND test_c_api_lp)

    add_executable(test_c_api_contract
        tests/unit/test_c_api_contract.cpp
    )
    target_link_libraries(test_c_api_contract PRIVATE processinglab_geometry)
    add_test(NAME c_api_contract COMMAND test_c_api_contract)

    add_executable(test_ub_from_psi
        tests/unit/test_ub_from_psi.cpp
    )
    target_link_libraries(test_ub_from_psi PRIVATE processinglab_geometry)
    add_test(NAME ub_from_psi COMMAND test_ub_from_psi)

    add_executable(test_lab_frame
        tests/unit/test_lab_frame.cpp
    )
    target_link_libraries(test_lab_frame PRIVATE processinglab_geometry)
    add_test(NAME lab_frame COMMAND test_lab_frame)

    if(MOSFLM_FORTRAN_TESTS)
        add_executable(test_fortran_c_api
            fortran/mosflm_c_bind.f90
            fortran/test_fortran_c_api.f90
        )
        target_link_libraries(test_fortran_c_api PRIVATE processinglab_geometry)
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
        target_link_libraries(test_predict_still_f77 PRIVATE processinglab_geometry)
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
        target_link_libraries(test_predict_cpp_bridge PRIVATE processinglab_geometry)
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
        target_link_libraries(test_predict_cpp_psi PRIVATE processinglab_geometry)
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
        target_link_libraries(test_predict_gendata_pack PRIVATE processinglab_geometry)
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
        target_link_libraries(test_predict_cpp_hook PRIVATE processinglab_geometry)
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
        target_link_libraries(test_predict_cpp_hook_f77 PRIVATE processinglab_geometry)
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
        target_link_libraries(test_predict_cpp_from_common PRIVATE processinglab_geometry)
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
        add_library(processinglab_predict_hook STATIC
            fortran/predict_cpp_bridge.f90
            fortran/predict_gendata_pack.f90
            fortran/predict_cpp_hook.f90
            fortran/predict_cpp_from_common.f90
        )
        target_link_libraries(processinglab_predict_hook PUBLIC processinglab_geometry)
        set_target_properties(processinglab_predict_hook PROPERTIES
            Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_hook_lib_mods
            POSITION_INDEPENDENT_CODE ON
        )
        message(STATUS "MOSFLM_FORTRAN_ORACLE=ON: building processinglab_predict_hook")
    elseif(MOSFLM_FORTRAN_ORACLE AND NOT MOSFLM_HAS_FORTRAN)
        message(WARNING
            "MOSFLM_FORTRAN_ORACLE=ON but no Fortran compiler; "
            "skipping processinglab_predict_hook")
    endif()

    add_executable(test_known_indexed_oracle
        tests/differential/test_known_indexed_oracle.cpp
    )
    target_link_libraries(test_known_indexed_oracle PRIVATE processinglab_geometry)
    add_test(
        NAME known_indexed_oracle
        COMMAND test_known_indexed_oracle
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_insulin_hkl_overlap
        tests/differential/test_insulin_hkl_overlap.cpp
    )
    target_link_libraries(test_insulin_hkl_overlap PRIVATE processinglab_geometry)
    add_test(
        NAME insulin_hkl_overlap
        COMMAND test_insulin_hkl_overlap
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_insulin_xdet
        tests/differential/test_insulin_xdet.cpp
    )
    target_link_libraries(test_insulin_xdet PRIVATE processinglab_geometry)
    add_test(
        NAME insulin_xdet
        COMMAND test_insulin_xdet
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_insulin_process_mtz
        tests/differential/test_insulin_process_mtz.cpp
    )
    target_link_libraries(test_insulin_process_mtz PRIVATE processinglab_geometry)
    add_test(
        NAME insulin_process_mtz
        COMMAND test_insulin_process_mtz
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_insulin_cpp_process
        tests/differential/test_insulin_cpp_process.cpp
    )
    target_link_libraries(test_insulin_cpp_process PRIVATE processinglab_geometry)
    add_test(
        NAME insulin_cpp_process
        COMMAND test_insulin_cpp_process
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_insulin_isi
        tests/differential/test_insulin_isi.cpp
    )
    target_link_libraries(test_insulin_isi PRIVATE processinglab_geometry)
    add_test(
        NAME insulin_isi
        COMMAND test_insulin_isi
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_insulin_profile
        tests/differential/test_insulin_profile.cpp
    )
    target_link_libraries(test_insulin_profile PRIVATE processinglab_geometry)
    add_test(
        NAME insulin_profile
        COMMAND test_insulin_profile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_insulin_postref
        tests/differential/test_insulin_postref.cpp
    )
    target_link_libraries(test_insulin_postref PRIVATE processinglab_geometry)
    add_test(
        NAME insulin_postref
        COMMAND test_insulin_postref
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_executable(test_insulin_control
        tests/differential/test_insulin_control.cpp
    )
    target_link_libraries(test_insulin_control PRIVATE processinglab_geometry)
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
        target_link_libraries(test_bravais_classifier PRIVATE processinglab_indexing)
        add_test(NAME bravais_classifier COMMAND test_bravais_classifier)

        add_executable(test_bravais_select
            tests/unit/test_bravais_select.cpp
        )
        target_link_libraries(test_bravais_select PRIVATE processinglab_indexing)
        add_test(
            NAME bravais_select
            COMMAND test_bravais_select
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_executable(test_spot_list
            tests/unit/test_spot_list.cpp
        )
        target_link_libraries(test_spot_list PRIVATE processinglab_indexing)
        add_test(
            NAME spot_list
            COMMAND test_spot_list
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_executable(test_spot_finder
            tests/unit/test_spot_finder.cpp
        )
        target_link_libraries(test_spot_finder PRIVATE processinglab_indexing)
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
        target_link_libraries(test_insulin_autoindex PRIVATE processinglab_indexing)
        add_test(
            NAME insulin_autoindex
            COMMAND test_insulin_autoindex
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        # Thaumatin Eiger HDF5 FINDSPOTS→AUTOINDEX (soft-skip without master).
        # Runs in the build dir: AUTOINDEX writes ./XPARM.XDS.
        add_executable(test_thaumatin_autoindex
            tests/differential/test_thaumatin_autoindex.cpp
        )
        target_link_libraries(test_thaumatin_autoindex PRIVATE processinglab_indexing)
        add_test(
            NAME thaumatin_autoindex
            COMMAND test_thaumatin_autoindex
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
        set_tests_properties(thaumatin_autoindex PROPERTIES TIMEOUT 600)

        add_executable(test_control_autoindex
            tests/differential/test_control_autoindex.cpp
        )
        target_link_libraries(test_control_autoindex PRIVATE processinglab_indexing)
        add_test(
            NAME control_autoindex
            COMMAND test_control_autoindex
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_executable(test_diff_vector_index
            tests/unit/test_diff_vector_index.cpp
        )
        target_link_libraries(test_diff_vector_index PRIVATE processinglab_indexing)
        add_test(
            NAME diff_vector_index
            COMMAND test_diff_vector_index
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_executable(test_supercell_reduce
            tests/unit/test_supercell_reduce.cpp
        )
        target_link_libraries(test_supercell_reduce PRIVATE processinglab_indexing)
        add_test(
            NAME supercell_reduce
            COMMAND test_supercell_reduce
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_executable(test_multiframe_spot_finder
            tests/unit/test_multiframe_spot_finder.cpp
        )
        target_link_libraries(test_multiframe_spot_finder PRIVATE processinglab_indexing)
        add_test(
            NAME multiframe_spot_finder
            COMMAND test_multiframe_spot_finder
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_executable(test_max_cell
            tests/unit/test_max_cell.cpp
        )
        target_link_libraries(test_max_cell PRIVATE processinglab_indexing)
        add_test(
            NAME max_cell
            COMMAND test_max_cell
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_executable(test_fft3d_index
            tests/unit/test_fft3d_index.cpp
        )
        target_link_libraries(test_fft3d_index PRIVATE processinglab_indexing)
        add_test(
            NAME fft3d_index
            COMMAND test_fft3d_index
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_executable(test_rlp_roundtrip
            tests/unit/test_rlp_roundtrip.cpp
        )
        target_link_libraries(test_rlp_roundtrip PRIVATE processinglab_indexing)
        add_test(
            NAME rlp_roundtrip
            COMMAND test_rlp_roundtrip
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_executable(test_idxref_refine
            tests/unit/test_idxref_refine.cpp
        )
        target_link_libraries(test_idxref_refine PRIVATE processinglab_indexing)
        add_test(
            NAME idxref_refine
            COMMAND test_idxref_refine
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_executable(test_76c1b_autoindex
            tests/differential/test_76c1b_autoindex.cpp
        )
        target_link_libraries(test_76c1b_autoindex PRIVATE processinglab_indexing)
        add_test(
            NAME 76c1b_autoindex
            COMMAND test_76c1b_autoindex
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
    endif()

    # 76c1b PROCESS I/σ + CORRECT-lite gates (XPARM-seeded; soft-skip if missing).
    add_executable(test_76c1b_process_quality
        tests/differential/test_76c1b_process_quality.cpp
    )
    target_link_libraries(test_76c1b_process_quality PRIVATE processinglab_geometry)
    add_test(
        NAME 76c1b_process_quality
        COMMAND test_76c1b_process_quality
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()
