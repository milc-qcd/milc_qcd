add_library(scidac INTERFACE)

set(HAVEQMP OFF)
set(HAVEQIO OFF)

set(QMP_TAG 3010fef)
set(QIO_TAG qio3-0-0)
if(MPP)
    set(QMP_OPTIONS "QMP_MPI ON")
    set(QIO_OPTIONS "QIO_ENABLE_PARALLEL_BUILD ON"
        "QIO_ENABLE_OUTPUT_BUFFERING ON" "QIO_ENABLE_FAST_ROUTE ON"
    )
else()
    set(QMP_OPTIONS "QMP_MPI OFF")
    set(QIO_OPTIONS "QIO_ENABLE_PARALLEL_BUILD OFF")
endif()

if(WANTQMP)
    if(DOWNLOAD_SCIDAC)
        CPMAddPackage(
            NAME QMP
            GITHUB_REPOSITORY usqcd-software/qmp
            GIT_TAG ${QMP_TAG}
            OPTIONS ${QMP_OPTIONS}
        )
    else()
        find_package(QMP REQUIRED)
    endif()
    set(HAVEQMP ON)
    if(MPP)
        target_compile_definitions(scidac INTERFACE QMP_MPI)
    else()
        target_compile_definitions(scidac INTERFACE QMP_SPI)
    endif()
    target_compile_definitions(scidac INTERFACE HAVE_QMP)
    target_link_libraries(scidac INTERFACE QMP::qmp)
endif()

if(WANTQIO)
    if(NOT WANTQMP)
        message(FATAL_ERROR "Use of QIO (via WANTQIO=ON) requires QMP. Please set WANTQMP=ON.")
    endif()
    if(DOWNLOAD_SCIDAC)
        CPMAddPackage(
            NAME QIO
            GITHUB_REPOSITORY usqcd-software/qio
            GIT_TAG ${QIO_TAG}
            OPTIONS ${QIO_OPTIONS}
        )
    else()
        find_package(QIO REQUIRED)
    endif()
    set(HAVEQIO ON)
    target_compile_definitions(scidac INTERFACE HAVE_QIO)
    target_link_libraries(scidac INTERFACE QIO::qio)
endif()
