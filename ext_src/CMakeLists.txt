# Source groups and targets follow Make_template.

add_interface_library(my_objects
    clover_info.c
    gauge_info.c
    ksprop_info.c
    ks_source_info.c
    make_ext_src.c
    setup.c
    w_source_info.c
)

add_interface_library(g_objects
    ../generic/layout_hyper_prime.c
)

add_interface_library(objects
    OBJECTS
    my_objects
    g_objects
    cl_io_objects
    fft_objects
    fermion_objects
    ks_io_objects
    standard_objects
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

add_milc_executable(ext_src control.c
    QUARK asqtad_action
    DEFINES FN HAVE_KS HAVE_DIRAC NO_GAUGE_FIELD
)

