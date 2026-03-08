add_interface_library(my_objects
    clover_info.c
    gauge_info.c
    ../generic_ks/io_prop_ks_fm.c
    ../generic_wilson/io_prop_w_fm.c
    ksprop_info.c
    ks_source_info.c
    w_source_info.c
    make_prop.c
    setup.c
    spectrum_cl.c
)

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
    cl_io_objects
    cl_objects
    cl_spectrum_objects
    fermion_objects
    fft_objects
    gauge_objects
    ks_io_objects
    layout
    standard_objects
    add_objects
)

add_interface_library(libraries
    OBJECTS
    qcd
    scidac
    add
    m
)

add_milc_executable(su3_clov control.c
    QUARK asqtad_action
    DEFINES HAVE_DIRAC FN HAVE_KS
    EXTRA_OBJECTS fn_objects
)

set(hisq_options
    FN
    HISQ_REUNIT_ALLOW_SVD
    HISQ_REUNIT_SVD_REL_ERROR=1e-8
    HISQ_REUNIT_SVD_ABS_ERROR=1e-8
)

add_milc_executable(su3_clov_hisq control.c
    QUARK hisq/hisq_u3_action
    DEFINES HAVE_DIRAC HAVE_KS ${hisq_options}
    EXTRA_OBJECTS hisq_objects
)
