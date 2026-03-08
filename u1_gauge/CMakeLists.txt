add_interface_library(my_objects
    setup.c
    gauge_info.c
)

add_interface_library(g_objects
    ../generic/field_utilities.c
    ../generic/io_detect.c
    ../generic/io_helpers.c
    ../generic/io_lat_utils.c
    ../generic/make_lattice.c
    ../generic/remap_stdio_from_args.c
    ../generic/ranstuff.c
    ../generic/remap_fftw_fourier.c
    ../generic/show_generic_opts.c
)

add_interface_library(g_u1_objects
    ../generic_u1/io_u1lat.c
    ../generic_u1/u1avlink.c
    ../generic_u1/u1link.c
    ../generic_u1/u1plaq.c
    ../generic_u1/u1pot.c
)

if(HAVEQIO)
    append_interface_library(g_objects
        ../generic/file_types_milc_usqcd.c
        ../generic/io_scidac.c
        ../generic/io_scidac_types.c
    )
endif()

add_interface_library(layout
    ../generic/layout_hyper_prime.c
)

add_interface_library(add_objects
    OBJECTS
    machine_dep_io
    commpkg
)

add_interface_library(objects
    OBJECTS
    my_objects
    g_objects
    g_u1_objects
    layout
    add_objects
)

add_interface_library(libraries
    OBJECTS
    qcd
    scidac
    add
    m
)

add_milc_executable(u1_g control.c
    DEFINES NO_GAUGE_FIELD GAUGE=COULOMB SCHEME=QED_L
    EXTRA_SOURCES momgauge.c
)

add_milc_executable(u1_g-Basak-Coul-QED_TL control.c
    DEFINES NO_GAUGE_FIELD GAUGE_METHOD=MILC_ORIGINAL
    EXTRA_SOURCES momgauge.c
)

add_milc_executable(u1_g-Hatton-Feynman-QED_L control.c
    DEFINES NO_GAUGE_FIELD GAUGE_METHOD=HATTON_REVISION
    EXTRA_SOURCES momgauge.c
)

add_milc_executable(u1_g-Feynman-QED_L control.c
    DEFINES NO_GAUGE_FIELD GAUGE=FEYNMAN SCHEME=QED_L
    EXTRA_SOURCES momgauge.c
)

add_milc_executable(u1_g-Lorenz-QED_L control.c
    DEFINES NO_GAUGE_FIELD GAUGE=LORENZ SCHEME=QED_L
    EXTRA_SOURCES momgauge.c
)

add_milc_executable(u1_g-Coulomb-QED_L control.c
    DEFINES NO_GAUGE_FIELD GAUGE=COULOMB SCHEME=QED_L
    EXTRA_SOURCES momgauge.c
)

add_milc_executable(u1_g-Coulomb-QED_TL control.c
    DEFINES NO_GAUGE_FIELD GAUGE=COULOMB SCHEME=QED_TL
    EXTRA_SOURCES momgauge.c
)

add_milc_executable(u1_g_convert control_convert.c
    DEFINES NO_GAUGE_FIELD
    EXTRA_SOURCES momgauge.c
)
