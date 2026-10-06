# Source groups and targets follow Make_template.

add_interface_library(my_objects
    accumulate_density.c
    gauge_info.c
    linearlsq.c
    print_corr.c
    rcorr.c
    setup.c
    symmetrize.c
)

add_interface_library(g_objects
    ../generic/layout_hyper_prime.c
)

add_interface_library(objects
    OBJECTS
    my_objects
    g_objects
    standard_objects
    fft_objects
    machine_dep_io
    commpkg
)

add_interface_library(libraries
    OBJECTS
    qcd
    scidac
    add
    m
)

add_milc_executable(current_rcorr control.c
    DEFINES NO_GAUGE_FIELD
)

