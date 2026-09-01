cmake_minimum_required(VERSION 3.18)
project(milc_qcd LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)

# Options
set(SUBDIR "ks_imp_rhmc" CACHE STRING "Project subdirectory")
set_property(CACHE SUBDIR PROPERTY STRINGS ks_imp_rhmc ks_spectrum)
set(ARCH "general" CACHE STRING "CPU architecture (knl, knc, skx, hsw, general)")
set_property(CACHE ARCH PROPERTY STRINGS knl knc skx hsw general)
set(OFFLOAD "cuda" CACHE STRING "Accelerator offload backend (cuda, hip, sycl, openmp)")
set_property(CACHE OFFLOAD PROPERTY STRINGS cuda hip sycl openmp)
option(MPP "Enable MPI" OFF)
option(OMP "Enable OpenMP" OFF)
set(PRECISION 2 CACHE STRING "Precision (1 = single, 2 = double)")
set_property(CACHE PRECISION PROPERTY STRINGS 1 2)
set(MIXED_PRECISION 2 CACHE STRING "Mixed-precision solver for GPU (1 = half, 2 = max)")
set_property(CACHE MIXED_PRECISION PROPERTY STRINGS 1 2)
set(MIXED_PRECISION_GPU ${MIXED_PRECISION} CACHE STRING "Mixed-precision for GPU (1 = half, 2 = max)")
set_property(CACHE MIXED_PRECISION_GPU PROPERTY STRINGS 1 2)
option(DOWNLOAD_SCIDAC "Download SciDAC libraries" OFF)
option(WANTQMP "Enable QMP communication library" OFF)
option(WANTQIO "Enable QIO binary I/O library (requires WANTQMP)" OFF)
option(WANT_APE_IO "Enable APE I/O support" OFF)
option(WANTFFTW "Enable FFTW" OFF)
option(WANTPRIMME "Enable PRIMME" OFF)
option(WANTARPACK "Enable ARPACK" OFF)
option(WANTQUDA "Enable QUDA" OFF)
option(WANTQPHIX "Enable QPhiX" OFF)
option(WANTHADRONS "Enable Hadrons (not yet wired in top-level CMake)" OFF)
option(WANTGRID "Enable Grid" OFF)
option(WANTQPHIXJ "Enable QPhiXJ (not yet wired in top-level CMake)" OFF)

option(GPU_DEFAULT "Default value for all GPU options" OFF)
mark_as_advanced(GPU_DEFAULT)
option(GPU_CL_BCG "Enable Clover BCG GPU" ${GPU_DEFAULT})
option(GPU_FN_CG "Enable KS FN CG GPU" ${GPU_DEFAULT})
option(GPU_GA "Enable gauge action GPU" ${GPU_DEFAULT})
option(GPU_GF "Enable gauge force GPU" ${GPU_DEFAULT})
option(GPU_FL "Enable fermion links GPU" ${GPU_DEFAULT})
option(GPU_FF "Enable fermion force GPU" ${GPU_DEFAULT})
option(GPU_EIG "Enable eigensolver GPU" ${GPU_DEFAULT})
option(GPU_GSMEAR "Enable gaussian smearing GPU" ${GPU_DEFAULT})
option(GPU_KS_CONT "Enable KS contraction GPU" ${GPU_DEFAULT})
option(GPU_SHIFT "Enable shift GPU" ${GPU_DEFAULT})
option(GPU_SPIN_TASTE "Enable spin taste GPU" ${GPU_DEFAULT})
option(GPU_GAUGEFIX_OVR "Enable gauge fixing override GPU" ${GPU_DEFAULT})
option(GPU_MULTIGRID "Enable multigrid GPU" ${GPU_DEFAULT})

option(QPHIX_DEFAULT "Default value for all QPhiX options" OFF)
mark_as_advanced(QPHIX_DEFAULT)
option(QPHIX_FN_CG "Enable KS FN CG QPhiX" ${QPHIX_DEFAULT})
option(QPHIX_GF "Enable gauge force QPhiX" ${QPHIX_DEFAULT})

option(GAUSS_SMEAR_KS_TWOLINK "Compute two-link field for KS Gaussian smearing" ON)

option(INLINE_C "Enable C inline functions" OFF)
option(INLINE_C_GLOBAL "Enable C global inline functions" ON)
option(INLINE_SSE "Enable SSE inline functions" OFF)
option(INLINE_SSE_GLOBAL "Enable SSE global inline functions" OFF)

option(TIME_DEFAULT "Default value for all timing options" ON)
mark_as_advanced(TIME_DEFAULT)
option(TIME_NERSC "Enable NERSC file timing" ${TIME_DEFAULT})
option(TIME_CG "Enable CG timing" ${TIME_DEFAULT})
option(TIME_FF "Enable fermion force timing" ${TIME_DEFAULT})
option(TIME_FL "Enable fermion links timing" ${TIME_DEFAULT})
option(TIME_GA "Enable gauge action timing" ${TIME_DEFAULT})
option(TIME_GF "Enable gauge force timing" ${TIME_DEFAULT})
option(TIME_REMAP "Enable remap timing" ${TIME_DEFAULT})
option(TIME_PR "Enable propagator timing" ${TIME_DEFAULT})
option(TIME_IO "Enable I/O timing" ${TIME_DEFAULT})
option(TIME_WM "Enable Wilson matrix timing" ${TIME_DEFAULT})
option(TIME_GS "Enable Gaussian smearing timing" ${TIME_DEFAULT})

option(DEBUG_CG_OK "Print inverter convergence information even when successful" ON)
option(DEBUG_REMAP_STDIO_APPEND "Append stdout from all MPI ranks for remapped stdio" ON)
option(DEBUG_CHECK_MALLOC "Report malloc/free activity for debugging" OFF)

option(COMPAT_QOPQDP_NORM "Use old QOPQDP norm" OFF)
option(COMPAT_STAGGERED2NAIVE "Use pre-7.7.2 staggered-to-naive conversion" OFF)
option(COMPAT_GAUSSRAND "Use old Gaussian random number generator" OFF)

option(GEOM_FIX_NODE "Fix node geometry" OFF)
option(GEOM_FIX_IONODE "Fix ionode geometry" OFF)

option(SHIFT_ONE_SIDED "Use forward-only shift for nonlocal KS operators" OFF)

option(CLMEM_LEAN "Write intermediate clover propagators to disk to save memory" OFF)

set(KS_MULTICG "HYBRID" CACHE STRING "KS multi-CG algorithm (OFFSET, HYBRID, FAKE, REVERSE, REVHYB)")
set_property(CACHE KS_MULTICG PROPERTY STRINGS OFFSET HYBRID FAKE REVERSE REVHYB)
set(KS_MULTIFF "FNMAT" CACHE STRING "KS multi-fermion force algorithm (FNMAT, FNMATREV, ASVEC)")
set_property(CACHE KS_MULTIFF PROPERTY STRINGS FNMAT FNMATREV ASVEC)
set(CL_CG "BICG" CACHE STRING "Clover inverter algorithm (BICG, CG, MR, HOP)")
set_property(CACHE CL_CG PROPERTY STRINGS BICG CG MR HOP)

include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/libraries/include
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(cmake/CPM.cmake)

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
    elseif(ARCH STREQUAL "icx")
        add_compile_options(-xCORE-AVX512)
        set(binext ".icx")
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
    set(HAVE_PRIMME ON)
    target_compile_definitions(add INTERFACE HAVE_PRIMME)
    target_link_libraries(add INTERFACE PRIMME::primme)
    # target_compile_definitions(add INTERFACE PRIMME_PRECOND MATVEC_PRECOND POLY_EIGEN)
endif()

if(WANTARPACK)
    find_package(ARPACK REQUIRED)
    set(HAVE_ARPACK ON)
    target_compile_definitions(add INTERFACE HAVE_ARPACK)
    target_link_libraries(add INTERFACE ARPACK::arpack)
endif()

if(WANTQUDA)
    find_package(QUDA REQUIRED)
    set(HAVE_QUDA ON)
    set(HAVE_GPU ON)
    target_link_libraries(add INTERFACE QUDA::quda)
    target_compile_definitions(add INTERFACE HAVE_QUDA SET_QUDA_SUMMARIZE)

    if(GPU_FN_CG)
        set(GPU_EIG ON)
    endif()

    if(GPU_CL_BCG)
        set(HAVE_CL_GPU ON)
        target_compile_definitions(add INTERFACE USE_CL_GPU)
    endif()
    if(GPU_FN_CG)
        set(HAVE_FN_CG_GPU ON)
        target_compile_definitions(add INTERFACE USE_CG_GPU)
    endif()
    if(GPU_GA)
        set(HAVE_GA_GPU ON)
        target_compile_definitions(add INTERFACE USE_GA_GPU)
    endif()
    if(GPU_GF)
        set(HAVE_GF_GPU ON)
        target_compile_definitions(add INTERFACE USE_GF_GPU)
    endif()
    if(GPU_FL)
        set(HAVE_FL_GPU ON)
        target_compile_definitions(add INTERFACE USE_FL_GPU)
    endif()
    if(GPU_FF)
        set(HAVE_FF_GPU ON)
        target_compile_definitions(add INTERFACE USE_FF_GPU)
    endif()
    if(GPU_EIG)
        set(HAVE_EIG_GPU ON)
        target_compile_definitions(add INTERFACE USE_EIG_GPU)
    endif()
    if(GPU_GSMEAR)
        set(HAVE_GSMEAR_GPU ON)
        target_compile_definitions(add INTERFACE USE_GSMEAR_GPU)
    endif()
    if(GPU_KS_CONT)
        set(HAVE_KS_CONT_GPU ON)
        target_compile_definitions(add INTERFACE USE_KS_CONT_GPU)
    endif()
    if(GPU_SHIFT)
        set(HAVE_SHIFT_GPU ON)
        target_compile_definitions(add INTERFACE USE_SHIFT_GPU)
    endif()
    if(GPU_SPIN_TASTE)
        set(HAVE_SPIN_TASTE_GPU ON)
        target_compile_definitions(add INTERFACE USE_SPIN_TASTE_GPU)
    endif()
    if(GPU_GAUGEFIX_OVR)
        set(HAVE_GAUGEFIX_OVR_GPU ON)
        target_compile_definitions(add INTERFACE USE_GAUGEFIX_OVR_GPU)
    endif()
    if(GPU_MULTIGRID)
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

    if(QPHIX_FN_CG)
        set(HAVE_FN_CG_QPHIX ON)
        target_compile_definitions(add INTERFACE USE_CG_QPHIX)
    endif()
    if(QPHIX_GF)
        set(HAVE_GF_QPHIX ON)
        target_compile_definitions(add INTERFACE USE_GF_QPHIX)
    endif()
endif()

# if(WANTHADRONS)
#     set(WANTGRID ON)
#     set(HAVE_HADRONS ON)
#     target_compile_definitions(add INTERFACE HAVE_HADRONS)
# endif()

if(WANTGRID)
    find_package(Grid REQUIRED)
    set(HAVE_GRID ON)
    set(HAVE_GPU ON)
    target_link_libraries(add INTERFACE Grid::Grid)
    target_compile_definitions(add INTERFACE HAVE_GRID)

    if(GPU_FN_CG)
        set(HAVE_FN_CG_GPU ON)
        target_compile_definitions(add INTERFACE USE_CG_GPU)
    endif()
    if(GPU_FL)
        set(HAVE_FL_GPU ON)
        target_compile_definitions(add INTERFACE USE_FL_GPU)
    endif()
    if(GPU_FF)
        set(HAVE_FF_GPU ON)
        target_compile_definitions(add INTERFACE USE_FF_GPU)
    endif()
    if(GPU_GF)
        set(HAVE_GF_GPU ON)
        target_compile_definitions(add INTERFACE USE_GF_GPU)
    endif()
    if(GPU_EIG)
        set(HAVE_EIG_GPU ON)
        target_compile_definitions(add INTERFACE USE_EIG_GPU)
    endif()

    set(GRID_SHMEM_MAX 2048) # Megabytes
    set(GRID_DEVICE_MEM_MAX 32768) # Megabytes
    set(GRID_SHMEM_MPI 1)
    set(GRID_ACCELERATOR_THREADS 8)
    set(GRID_MULTI_CG GRID_5DCG) # GRID_5DCG GRID_BLOCKCG GRID_MRHSCG

    target_compile_definitions(add INTERFACE
        GRID_SHMEM_MAX=${GRID_SHMEM_MAX}
        GRID_DEVICE_MEM_MAX=${GRID_DEVICE_MEM_MAX}
        GRID_SHMEM_MPI=${GRID_SHMEM_MPI}
        GRID_ACCELERATOR_THREADS=${GRID_ACCELERATOR_THREADS}
        GRID_COMMS_OVERLAP
        GRID_MULTI_CG=${GRID_MULTI_CG}
    )
endif()

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
    $<$<BOOL:${INLINE_C}>:INLINE>
    $<$<BOOL:${INLINE_C_GLOBAL}>:C_GLOBAL_INLINE>
    $<$<BOOL:${INLINE_SSE}>:SSE_INLINE>
    $<$<BOOL:${INLINE_SSE_GLOBAL}>:SSE_GLOBAL_INLINE>
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
    $<$<BOOL:${TIME_NERSC}>:NERSC_TIME>
    $<$<BOOL:${TIME_CG}>:CGTIME>
    $<$<BOOL:${TIME_FF}>:FFTIME>
    $<$<BOOL:${TIME_FL}>:FLTIME>
    $<$<BOOL:${TIME_GA}>:GATIME>
    $<$<BOOL:${TIME_GF}>:GFTIME>
    $<$<BOOL:${TIME_REMAP}>:REMAP>
    $<$<BOOL:${TIME_PR}>:PRTIME>
    $<$<BOOL:${TIME_IO}>:IOTIME>
    $<$<BOOL:${TIME_WM}>:WMTIME>
    $<$<BOOL:${TIME_GS}>:GS_TIME>
)

set(cdebug
    $<$<BOOL:${DEBUG_CG_OK}>:CG_OK>
    $<$<BOOL:${DEBUG_REMAP_STDIO_APPEND}>:REMAP_STDIO_APPEND>
    $<$<BOOL:${DEBUG_CHECK_MALLOC}>:CHECK_MALLOC>
)

set(ccompat
    $<$<BOOL:${COMPAT_QOPQDP_NORM}>:OLD_QOPQDP_NORM>
    $<$<BOOL:${COMPAT_STAGGERED2NAIVE}>:OLD_STAGGERED2NAIVE>
    $<$<BOOL:${COMPAT_GAUSSRAND}>:OLD_GAUSSRAND>
)

set(cgeom
    $<$<BOOL:${GEOM_FIX_NODE}>:FIX_NODE_GEOM>
    $<$<BOOL:${GEOM_FIX_IONODE}>:FIX_IONODE_GEOM>
)

if(WANTQUDA)
    set(kscgstore FEWSUMS)
else()
    set(kscgstore DBLSTORE_FN FEWSUMS D_FN_GATHER13)
endif()

# set(kscgmulti
#     KS_MULTICG=${KS_MULTICG}
# )

# set(ksffmulti
#     KS_MULTIFF=${KS_MULTIFF}
# )

set(ksshift
    $<$<BOOL:${SHIFT_ONE_SIDED}>:ONE_SIDED_SHIFT>
)

# set(clcg
#     CL_CG=${CL_CG}
# )

set(clmem
    $<$<BOOL:${CLMEM_LEAN}>:CLOV_LEAN>
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
    ${ctime}
    ${inlineopt}
    # ${cprof}
    ${cdebug}
    ${ccompat}
    ${cgeom}
    ${kscgstore}
    # ${cprefetch}
    KS_MULTICG=${KS_MULTICG}
    # ${kscgmixed}
    # ${kssource}
    KS_MULTIFF=${KS_MULTIFF}
    ${ksshift}
    CL_CG=${CL_CG}
    ${clmem}
)

if(${MIXED_PRECISION_GPU} STREQUAL 1)
    target_compile_definitions(add INTERFACE HALF_MIXED)
elseif(${MIXED_PRECISION_GPU} STREQUAL 2)
    target_compile_definitions(add INTERFACE MAX_MIXED)
endif()

function(add_milc_executable TARGET MAIN)
    add_executable(${TARGET} ${MAIN})
    cmake_parse_arguments(AME "" "GAUGE;QUARK" "DEFINES;EXTRA_SOURCES;EXTRA_OBJECTS" ${ARGN})
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
    target_sources(${TARGET} PRIVATE ${AME_EXTRA_SOURCES} ${AME_UNPARSED_ARGUMENTS})
    target_link_libraries(${TARGET} PRIVATE objects ${AME_EXTRA_OBJECTS} libraries)
endfunction()

add_subdirectory(libraries)
add_subdirectory(${SUBDIR})
