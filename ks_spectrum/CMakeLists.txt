add_interface_library(my_objects
    SOURCES
    gauge_info.c
    ks_source_info.c
    ksprop_info.c
    make_prop.c
    setup.c
    spectrum_ks.c
)


add_interface_library(g_u1_objects
    SOURCES
    ../generic_u1/io_u1lat.c
    ../generic_u1/u1link.c
    ../generic_u1/u1plaq.c
)

add_interface_library(layout
    SOURCES
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
    gauge_objects
    ks_io_objects
    layout
    fft_objects
    fermion_objects
    ks_spectrum_objects
    standard_objects
    eigen_objects
    add_objects
)

add_interface_library(eigcg_objects
    SOURCES
    ../generic_ks/inc_eigcg.c
)

add_interface_library(libraries
    OBJECTS
    qcd
    scidac
    add
    m
    crypto
)

set(hisq_options
    FN
    HAVE_KS
    HISQ_REUNIT_ALLOW_SVD
    HISQ_REUNIT_SVD_REL_ERROR=1e-8
    DHISQ_REUNIT_SVD_ABS_ERROR=1e-8
)

add_milc_executable(ks_spectrum_hisq control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options}
    EXTRA_OBJECTS hisq_objects
)
