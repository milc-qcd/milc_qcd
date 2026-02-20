cmake_minimum_required(VERSION 3.18)
project(milc_qcd LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)

# Options
set(SUBDIR "ks_imp_rhmc" CACHE STRING "Project subdirectory")
set_property(CACHE SUBDIR PROPERTY STRINGS ks_imp_rhmc ks_spectrum)
set(ARCH "general" CACHE STRING "CPU architecture (knl, knc, skx, hsw, general)")
set_property(CACHE ARCH PROPERTY STRINGS knl knc skx hsw general)
set(OFFLOAD "cuda" CACHE STRING "GPU compiler (cuda, hip, sycl)")
set_property(CACHE OFFLOAD PROPERTY STRINGS cuda hip sycl)
option(MPP "Enable MPI" OFF)
option(OMP "Enable OpenMP" OFF)
set(PRECISION 2 CACHE STRING "Precision (1 = single, 2 = double)")
set_property(CACHE PRECISION PROPERTY STRINGS 1 2)
set(WANT_MIXED_PRECISION 2 CACHE STRING "Mixed-precision (1 = half, 2 = max)")
set_property(CACHE WANT_MIXED_PRECISION PROPERTY STRINGS 1 2)
set(WANT_MIXED_PRECISION_GPU ${WANT_MIXED_PRECISION} CACHE STRING "Mixed-precision for GPU (1 = half, 2 = max)")
set_property(CACHE WANT_MIXED_PRECISION_GPU PROPERTY STRINGS 1 2)
option(WANTQOP "Enable QOP" OFF)
option(WANTQIO "Enable QIO" OFF)
option(WANTQMP "Enable QMP" OFF)
option(WANT_APE_IO "Enable APE IO" OFF)
option(WANTFFTW "Enable FFTW" OFF)
option(WANTPRIMME "Enable PRIMME" OFF)
option(WANTARPACK "Enable ARPACK" OFF)
option(WANTQUDA "Enable QUDA" OFF)
option(WANTQPHIX "Enable QPhiX" OFF)
option(WANTHADRONS "Enable Hadrons" OFF)
option(WANTGRID "Enable Grid" OFF)
option(WANTQPHIXJ "Enable QPhiXJ" OFF)

option(WANT_CL_BCG_GPU "Enable Clover BCG GPU" OFF)
option(WANT_FN_CG_GPU "Enable KS FN CG GPU" OFF)
option(WANT_GA_GPU "Enable gauge action CG GPU" OFF)
option(WANT_GF_GPU "Enable gauge force GPU" OFF)
option(WANT_FL_GPU "Enable fermion links GPU" OFF)
option(WANT_FF_GPU "Enable fermion force GPU" OFF)
option(WANT_EIG_GPU "Enable eigensolver GPU" OFF)
option(WANT_GSMEAR_GPU "Enable gaussian smearing GPU" OFF)
option(WANT_KS_CONT_GPU "Enable KS contraction GPU" OFF)
option(WANT_SHIFT_GPU "Enable shift GPU" OFF)
option(WANT_SPIN_TASTE_GPU "Enable spin taste GPU" OFF)
option(WANT_GAUGEFIX_OVR_GPU "Enable gauge fixing override GPU" OFF)
option(WANT_MULTIGRID "Enable multigrid GPU" OFF)

option(WANT_FN_CG_QPHIX "Enable KS FN CG QPhiX" OFF)
option(WANT_GF_QPHIX "Enable gauge force QPhiX" OFF)

option(GAUSS_SMEAR_KS_TWOLINK "Use two-link term in Gaussian smearing for KS" OFF)

include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/libraries/include
)

function(add_interface_library TARGET)
    cmake_parse_arguments(AIL "" "" "OBJECTS;SOURCES" ${ARGN})
    add_library(${TARGET} INTERFACE)
    target_sources(${TARGET} INTERFACE ${AIL_UNPARSED_ARGUMENTS} ${AIL_SOURCES})
    target_link_libraries(${TARGET} INTERFACE ${AIL_OBJECTS})
endfunction()

function(append_interface_library TARGET)
    cmake_parse_arguments(AIL "" "" "OBJECTS;SOURCES" ${ARGN})
    target_sources(${TARGET} INTERFACE ${AIL_UNPARSED_ARGUMENTS} ${AIL_SOURCES})
    target_link_libraries(${TARGET} INTERFACE ${AIL_OBJECTS})
endfunction()

# Intel Classic Compiler Configuration
if(CMAKE_C_COMPILER_ID STREQUAL "Intel" OR CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")
    if(ARCH STREQUAL "knl")
        add_compile_options(-xMIC-AVX512)
        set(binext ".knl")
    elseif(ARCH STREQUAL "knc")
        add_compile_options(-mmic)
        set(binext ".knc")
    elseif(ARCH STREQUAL "skx")
        add_compile_options(-xCORE-AVX512 -qopt-zmm-usage=high)
        set(binext ".skx")
    elseif(ARCH STREQUAL "hsw")
        add_compile_options(-xCORE-AVX2)
        set(binext ".hsw")
    else()
        add_compile_options(-mavx)
        set(binext "")
    endif()

    if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
        add_compile_options(-parallel-source-info=2 -debug inline-debug-info -qopt-report=5)
    elseif(CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")
        add_compile_options(-parallel-source-info=2 -debug inline-debug-info -fsave-optimization-record)
    endif()
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "NVHPC")
    add_compile_options(-tp host)
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "IBM")
    add_compile_definitions(_LARGE_FILES)
else()
    add_compile_definitions(_FILE_OFFSET_BITS=64 _LARGEFILE64_SOURCE)
endif()

add_library(machine_dep_io INTERFACE)
target_sources(machine_dep_io INTERFACE
    generic/io_ansi.c
    # generic/io_nonansi.c
    # generic/io_dcap.c
)

if(MPP)
    add_compile_definitions(HAVE_MPI)
endif()

include(template_scidac.cmake)

add_library(add INTERFACE)

if(CMAKE_C_COMPILER_ID STREQUAL "Intel" OR CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")
    find_package(MKL REQUIRED)
    target_link_libraries(add INTERFACE MKL::MKL)
endif()

if(WANTFFTW)
    find_package(FFTW3 REQUIRED)
    set(HAVEFFTW ON)
    if(PRECISION STREQUAL 1)
        target_link_libraries(add INTERFACE FFTW3::fftw3f)
    elseif(PRECISION STREQUAL 2)
        target_link_libraries(add INTERFACE FFTW3::fftw3)
    endif()
endif()

if(WANTPRIMME)
    find_package(PRIMME REQUIRED)
    set(HAVEPRIMME ON)
    target_compile_definitions(add INTERFACE HAVE_PRIMME)
    target_link_libraries(add INTERFACE PRIMME::primme)
    # target_compile_definitions(add INTERFACE PRIMME_PRECOND POLY_EIGEN)
endif()

if(WANTARPACK)
    find_package(ARPACK REQUIRED)
    set(HAVEARPACK ON)
    target_compile_definitions(add INTERFACE HAVE_ARPACK)
    target_link_libraries(add INTERFACE ARPACK::arpack)
endif()

if(WANTQUDA)
    find_package(QUDA REQUIRED)
    set(HAVE_QUDA ON)
    set(HAVE_GPU ON)
    target_link_libraries(add INTERFACE QUDA::quda)
    target_compile_definitions(add INTERFACE HAVE_QUDA SET_QUDA_SUMMARIZE)

    if(WANT_FN_CG_GPU)
        set(WANT_EIG_GPU ON)
    endif()

    if(WANT_CL_BCG_GPU)
        set(HAVE_CL_GPU ON)
        target_compile_definitions(add INTERFACE USE_CL_GPU)
    endif()
    if(WANT_FN_CG_GPU)
        set(HAVE_FN_CG_GPU ON)
        target_compile_definitions(add INTERFACE USE_CG_GPU)
    endif()
    if(WANT_GA_GPU)
        set(HAVE_GA_GPU ON)
        target_compile_definitions(add INTERFACE USE_GA_GPU)
    endif()
    if(WANT_GF_GPU)
        set(HAVE_GF_GPU ON)
        target_compile_definitions(add INTERFACE USE_GF_GPU)
    endif()
    if(WANT_FL_GPU)
        set(HAVE_FL_GPU ON)
        target_compile_definitions(add INTERFACE USE_FL_GPU)
    endif()
    if(WANT_FF_GPU)
        set(HAVE_FF_GPU ON)
        target_compile_definitions(add INTERFACE USE_FF_GPU)
    endif()
    if(WANT_EIG_GPU)
        set(HAVE_EIG_GPU ON)
        target_compile_definitions(add INTERFACE USE_EIG_GPU)
    endif()
    if(WANT_GSMEAR_GPU)
        set(HAVE_GSMEAR_GPU ON)
        target_compile_definitions(add INTERFACE USE_GSMEAR_GPU)
    endif()
    if(WANT_KS_CONT_GPU)
        set(HAVE_KS_CONT_GPU ON)
        target_compile_definitions(add INTERFACE USE_KS_CONT_GPU)
    endif()
    if(WANT_SHIFT_GPU)
        set(HAVE_SHIFT_GPU ON)
        target_compile_definitions(add INTERFACE USE_SHIFT_GPU)
    endif()
    if(WANT_SPIN_TASTE_GPU)
        set(HAVE_SPIN_TASTE_GPU ON)
        target_compile_definitions(add INTERFACE USE_SPIN_TASTE_GPU)
    endif()
    if(WANT_GAUGEFIX_OVR_GPU)
        set(HAVE_GAUGEFIX_OVR_GPU ON)
        target_compile_definitions(add INTERFACE USE_GAUGEFIX_OVR_GPU)
    endif()
    if(WANT_MULTIGRID)
        target_compile_definitions(add INTERFACE MULTIGRID)
    endif()
endif()

if(WANTQPHIX)
    find_package(QPHIX REQUIRED)
    set(HAVE_QPHIX ON)
    target_compile_definitions(add INTERFACE HAVE_QPHIX QPHIX_PrecisionInt=${PRECISION})
    if(MPP)
        if(${ARCH} STREQUAL knl)
            target_link_libraries(add INTERFACE QPHIX::qphixmilc_avx512)
        elseif(${ARCH} STREQUAL knc)
            target_link_libraries(add INTERFACE QPHIX::qphixmilc_mic)
        elseif(${ARCH} STREQUAL hsw)
            target_link_libraries(add INTERFACE QPHIX::qphixmilc_avx2)
        elseif(${ARCH} STREQUAL skx)
            target_link_libraries(add INTERFACE QPHIX::qphixmilc_skx)
        endif()
    else()
        if(${ARCH} STREQUAL knl)
            target_link_libraries(add INTERFACE QPHIX::qphixmilc_avx512_single)
        elseif(${ARCH} STREQUAL knc)
            target_link_libraries(add INTERFACE QPHIX::qphixmilc_mic_single)
        elseif(${ARCH} STREQUAL hsw)
            target_link_libraries(add INTERFACE QPHIX::qphixmilc_avx2_single)
        elseif(${ARCH} STREQUAL skx)
            target_link_libraries(add INTERFACE QPHIX::qphixmilc_skx_single)
        endif()
    endif()

    if(WANT_FN_CG_QPHIX)
        set(HAVE_FN_CG_QPHIX ON)
        target_compile_definitions(add INTERFACE USE_CG_QPHIX)
    endif()
    if(WANT_GF_QPHIX)
        set(HAVE_GF_QPHIX ON)
        target_compile_definitions(add INTERFACE USE_GF_QPHIX)
    endif()
endif()

# if(WANTHADRONS)
#     set(WANTGRID ON)
#     set(HAVE_HADRONS ON)
#     target_compile_definitions(add INTERFACE HAVE_HADRONS)
# endif()

# if(WANTGRID)
#     find_package(Grid REQUIRED)
#     set(HAVE_GRID ON)
#     set(HAVE_GPU ON)
#     target_compile_definitions(add INTERFACE HAVE_GRID)
# endif()

# if(WANTQPHIXJ)
#     set(HAVE_QPHIXJ ON)
#     target_compile_definitions(add INTERFACE HAVE_QPHIXJ)
# endif()

if(WANT_APE_IO)
    set(HAVE_APE_IO ON)
    target_compile_definitions(add INTERFACE APE_LINKS_FILE)
endif()

if(GAUSS_SMEAR_KS_TWOLINK)
    target_compile_definitions(add INTERFACE GAUSS_SMEAR_KS_TWOLINK)
endif()

include(template_combos.cmake)

set(inlineopt
    # C_INLINE
    C_GLOBAL_INLINE
    # SSE_INLINE
    # SSE_GLOBAL_INLINE
)

find_package(Git REQUIRED)
execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --abbrev=4 --dirty --always --tags
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GITVERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
set(cgitver
    MILC_CODE_VERSION=\"${GITVERSION}\"
)

set(ctime
    NERSC_TIME CGTIME FFTIME FLTIME GFTIME REMAP PRTIME IOTIME WMTIME
)

set(cprof
)

set(cdebug
    CG_OK REMAP_STDIO_APPEND
)

set(ccompat)
# list(APPEND ccompat OLD_QOPQDP_NORM)
# list(APPEND ccompat OLD_STAGGERED2NAIVE OLD_GAUSSRAND)

set(cgeom
    # FIX_NODE_GEOM
    # FIX_IONODE_GEOM
)

if(WANTQUDA)
    set(kscgstore FEWSUMS)
else()
    set(kscgstore DBLSTORE_FN FEWSUMS D_FN_GATHER13)
endif()

set(cprefetch
)

set(kscgmulti
    KS_MULTICG=HYBRID
)

set(kscgmixed
)

set(kssource
)

set(ksffmulti
    KS_MULTIFF=FNMAT
    # KS_MULTIFF=FNMATREV
    # KS_MULTIFF=ASVEC
)

if(OMP)
    find_package(OpenMP REQUIRED)
    target_compile_definitions(add INTERFACE OMP)
    target_link_libraries(add INTERFACE OpenMP::OpenMP_C)
endif()

add_library(commpkg INTERFACE)
if(MPP)
    if(HAVEQMP)
        target_sources(commpkg INTERFACE generic/com_qmp.c)
        target_compile_definitions(commpkg INTERFACE QMP)
        target_link_libraries(commpkg INTERFACE QMP::qmp)
    else()
        find_package(MPI REQUIRED)
        target_sources(commpkg INTERFACE generic/com_mpi.c)
        target_compile_definitions(commpkg INTERFACE MPI_COMMS)
        target_link_libraries(commpkg INTERFACE MPI::MPI_C)
    endif()
else()
    if(HAVEQMP)
        target_sources(commpkg INTERFACE generic/com_qmp.c)
        target_compile_definitions(commpkg INTERFACE SINGLE)
        target_link_libraries(commpkg INTERFACE QMP::qmp)
    else()
        target_sources(commpkg INTERFACE generic/com_vanilla.c)
        target_compile_definitions(commpkg INTERFACE SINGLE)
    endif()
endif()

target_compile_definitions(add INTERFACE
    MILC_PRECISION=${PRECISION}
    ${cgitver}
    ${cscidac}
    ${ctime}
    ${cprof}
    ${cdebug}
    ${cgeom}
    ${kscgstore}
    ${cprefetch}
    ${kscgmulti}
    ${kscgmixed}
    ${kssource}
    ${ksffmulti}
    ${ccompat}
    ${cgitver}
)

if(${WANT_MIXED_PRECISION_GPU} STREQUAL 1)
    target_compile_definitions(add INTERFACE HALF_MIXED)
elseif(${WANT_MIXED_PRECISION_GPU} STREQUAL 2)
    target_compile_definitions(add INTERFACE MAX_MIXED)
endif()

function(add_milc_executable TARGET MAIN)
    add_executable(${TARGET} ${MAIN})
    cmake_parse_arguments(AME "" "GAUGE;QUARK" "DEFINES;EXTRA_OBJECTS" ${ARGN})
    set(gauge_include ${CMAKE_SOURCE_DIR}/generic/imp_actions)
    set(quark_include ${CMAKE_SOURCE_DIR}/generic_ks/imp_actions)
    set(local_include ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}_include)

    if(AME_GAUGE)
        set(SRC ${gauge_include}/${AME_GAUGE}.h)
        set(DST ${local_include}/gauge_action.h)
        add_custom_command(
            OUTPUT ${DST}
            COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC} ${DST}
            DEPENDS ${SRC}
            COMMENT "Copying ${AME_GAUGE}.h to gauge_action.h"
        )
        target_sources(${TARGET} PRIVATE ${DST})
    endif()

    if(AME_QUARK)
        set(SRC ${quark_include}/${AME_QUARK}.h)
        set(DST ${local_include}/quark_action.h)
        add_custom_command(
            OUTPUT ${DST}
            COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC} ${DST}
            DEPENDS ${SRC}
            COMMENT "Copying ${AME_QUARK}.h to quark_action.h"
        )
        target_sources(${TARGET} PRIVATE ${DST})
    endif()

    target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${local_include})
    target_compile_definitions(${TARGET} PRIVATE ${AME_DEFINES})
    target_link_libraries(${TARGET} PRIVATE objects ${AME_EXTRA_OBJECTS} libraries)
endfunction()

add_subdirectory(libraries)
add_subdirectory(${SUBDIR})
