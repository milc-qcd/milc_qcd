# Prefer a package-provided target. Older Grid installations use grid-config.
find_package(Grid CONFIG QUIET NO_MODULE)
if(TARGET Grid::Grid)
    set(Grid_FOUND TRUE)
    set(GRID_FOUND TRUE)
    return()
endif()

find_program(GRID_CONFIG_EXECUTABLE
    NAMES grid-config
    HINTS ${Grid_ROOT} ${GRID_ROOT} ${GRID_DIR} $ENV{Grid_ROOT} $ENV{GRID_ROOT} $ENV{GRID_DIR}
    PATH_SUFFIXES bin
    DOC "Path to the grid-config script"
)
mark_as_advanced(GRID_CONFIG_EXECUTABLE)

set(_grid_config_ok FALSE)
set(_grid_failure "grid-config was not found")
if(GRID_CONFIG_EXECUTABLE)
    set(_grid_config_ok TRUE)
    foreach(_query cxxflags ldflags libs prefix)
        execute_process(
            COMMAND "${GRID_CONFIG_EXECUTABLE}" --${_query}
            RESULT_VARIABLE _status
            OUTPUT_VARIABLE _grid_${_query}
            ERROR_VARIABLE _error
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(NOT _status STREQUAL "0")
            set(_grid_config_ok FALSE)
            set(_grid_failure "grid-config --${_query} failed: ${_error}")
            break()
        endif()
    endforeach()
endif()

# Empty compile/link flags are valid. A prefix is not a version number.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Grid
    REQUIRED_VARS GRID_CONFIG_EXECUTABLE _grid_config_ok _grid_libs
    REASON_FAILURE_MESSAGE "${_grid_failure}"
)
set(GRID_FOUND "${Grid_FOUND}")
if(NOT Grid_FOUND)
    return()
endif()
set(GRID_PREFIX "${_grid_prefix}")

# Keep the complete compiler flag string, including paired flags such as
# -isystem <path>, -include <header>, and accelerator options. SHELL preserves
# their grouping during CMake's option de-duplication. Grid flags are C++ only.
add_library(Grid::Grid INTERFACE IMPORTED)
if(_grid_cxxflags)
    target_compile_options(Grid::Grid INTERFACE
        "$<$<COMPILE_LANGUAGE:CXX>:SHELL:${_grid_cxxflags}>"
    )
endif()
if(_grid_ldflags)
    target_link_options(Grid::Grid INTERFACE "SHELL:${_grid_ldflags}")
endif()

# Libraries must follow object files at link time. Preserve their order and
# absolute archive paths instead of retaining only entries beginning with -l.
separate_arguments(_grid_libraries UNIX_COMMAND "${_grid_libs}")
set(_grid_link_libraries)
set(_grid_framework FALSE)
foreach(_item IN LISTS _grid_libraries)
    if(_grid_framework)
        list(APPEND _grid_link_libraries "-framework ${_item}")
        set(_grid_framework FALSE)
    elseif(_item STREQUAL "-framework")
        set(_grid_framework TRUE)
    else()
        list(APPEND _grid_link_libraries "${_item}")
    endif()
endforeach()
if(_grid_framework)
    message(FATAL_ERROR "grid-config --libs ends with -framework without a name")
endif()
target_link_libraries(Grid::Grid INTERFACE ${_grid_link_libraries})
