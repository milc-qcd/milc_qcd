set(GRID_ROOT_HINTS ${GRID_ROOT} ${GRID_DIR} $ENV{GRID_ROOT} $ENV{GRID_DIR})

find_program(GRID_CONFIG_EXECUTABLE
    NAMES grid-config
    HINTS ${GRID_ROOT_HINTS}
    PATH_SUFFIXES bin
    DOC "Path to the grid-config script"
)
mark_as_advanced(GRID_CONFIG_EXECUTABLE)

if(GRID_CONFIG_EXECUTABLE)
    execute_process(COMMAND ${GRID_CONFIG_EXECUTABLE} --cxxflags OUTPUT_VARIABLE _cxxflags OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${GRID_CONFIG_EXECUTABLE} --ldflags OUTPUT_VARIABLE _ldflags OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${GRID_CONFIG_EXECUTABLE} --libs OUTPUT_VARIABLE _libs OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${GRID_CONFIG_EXECUTABLE} --prefix OUTPUT_VARIABLE GRID_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)

    separate_arguments(_flags_list NATIVE_COMMAND "${_cxxflags}")
    foreach(_flag ${_flags_list})
        if(_flag MATCHES "^-D")
            string(REGEX REPLACE "^-D" "" _definitions "${_flag}")
            list(APPEND GRID_DEFINITIONS "${_definitions}")
        elseif(_flag MATCHES "^-I")
            string(REGEX REPLACE "^-I" "" _include_dirs "${_flag}")
            list(APPEND GRID_INCLUDE_DIRS "${_include_dirs}")
        endif()
    endforeach()

    separate_arguments(_flags_list NATIVE_COMMAND "${_ldflags}")
    foreach(_flag ${_flags_list})
        if(_flag MATCHES "^-L")
            string(REGEX REPLACE "^-L" "" _library_dirs "${_flag}")
            list(APPEND GRID_LIBRARY_DIRS "${_library_dirs}")
        endif()
    endforeach()

    separate_arguments(_flags_list NATIVE_COMMAND "${_libs}")
    foreach(_flag ${_flags_list})
        if(_flag MATCHES "^-l")
            string(REGEX REPLACE "^-l" "" _libraries "${_flag}")
            list(APPEND GRID_LIBRARIES "${_libraries}")
        endif()
    endforeach()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Grid
    REQUIRED_VARS GRID_CONFIG_EXECUTABLE GRID_DEFINITIONS GRID_INCLUDE_DIRS GRID_LIBRARY_DIRS GRID_LIBRARIES
    VERSION_VAR GRID_PREFIX
)

if(GRID_FOUND AND NOT TARGET Grid::Grid)
    add_library(Grid::Grid INTERFACE IMPORTED)
    set_target_properties(Grid::Grid PROPERTIES
        INTERFACE_COMPILE_DEFINITIONS "${GRID_DEFINITIONS}"
        INTERFACE_INCLUDE_DIRECTORIES "${GRID_INCLUDE_DIRS}"
        INTERFACE_LINK_DIRECTORIES "${GRID_LIBRARY_DIRS}"
        INTERFACE_LINK_LIBRARIES "${GRID_LIBRARIES}"
    )
endif()
