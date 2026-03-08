add_interface_library(my_objects
    setup.c
    gauge_info.c
    io_helpers_ks_utils.c
    ks_info.c
)

add_interface_library(g_objects
    # ../generic/gauge_stuff.c
    # ../generic/ploop3.c
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
    add_objects
    fermion_objects
    fft_objects
    g_objects
    gauge_objects
    ks_io_objects
    layout
    my_objects
    standard_objects
)

add_interface_library(libraries
    OBJECTS
    qcd
    scidac
    add
    m
)

add_milc_executable(check_ks_invert control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES FN HAVE_KS CHECK_INVERT
    EXTRA_OBJECTS fn_objects
    EXTRA_SOURCES check_invert2.c check_ks_invert.c
)

add_milc_executable(check_ks_invert_hisq control.c
    GAUGE symanzik_1loop_action
    QUARK hisq/hisq_u3_action
    DEFINES FN HAVE_KS CHECK_INVERT
    EXTRA_OBJECTS hisq_objects
    EXTRA_SOURCES check_invert2.c check_ks_invert.c
    # control_ks_invert.c setup_ks_invert.c
)

set(hisq_link_opts
    FN
    HAVE_KS
    HISQ_REUNIT_ALLOW_SVD
    HISQ_REUNIT_SVD_REL_ERROR=1e-8
    HISQ_REUNIT_SVD_ABS_ERROR=1e-8
)

add_milc_executable(make_links_hisq control.c
    GAUGE symanzik_1loop_action
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_link_opts} LINK_FATTENING
    EXTRA_OBJECTS hisq_objects hisq_force
    EXTRA_SOURCES check_link_fattening.c
)

add_milc_executable(check_link_fattening control.c
    GAUGE symanzik_1loop_action
    QUARK hisq/hisq_u3_test_action
    DEFINES ${hisq_link_opts} CHECK_FATTENING
    EXTRA_OBJECTS hisq_objects hisq_force
)

add_milc_executable(check_fermion_force_asqtad control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES FN HAVE_KS FERMION_FORCE
    EXTRA_OBJECTS fn_objects asq_force
    EXTRA_SOURCES
    check_fermion_force.c
    ../generic_ks/load_rhmc_params.c
    ../generic_ks/grsource_rhmc.c
    ../generic_ks/ks_ratinv.c
)

set(hisq_force_opts
    FN
    HAVE_KS
    FERMION_FORCE
    HISQ_REUNIT_ALLOW_SVD
    HISQ_REUNIT_SVD_REL_ERROR=1e-9
    HISQ_REUNIT_SVD_ABS_ERROR=1e-9
    HISQ_SVD_COUNTER
    HISQ_FF_MULTI_WRAPPER
    HISQ_FORCE_FILTER=5.0e-5
    HISQ_FORCE_FILTER_COUNTER
)

add_milc_executable(check_fermion_force_hisq control.c
    GAUGE symanzik_1loop_action
    QUARK hisq/hisq_u3_test_action
    DEFINES ${hisq_force_opts}
    EXTRA_OBJECTS hisq_objects hisq_force
    EXTRA_SOURCES
    check_fermion_force.c
    ../generic_ks/load_rhmc_params.c
    ../generic_ks/grsource_rhmc.c
    ../generic_ks/ks_ratinv.c
)

add_milc_executable(check_reunit_deriv control.c
    GAUGE symanzik_1loop_action
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_force_opts} REUNIT
    EXTRA_OBJECTS hisq_objects hisq_force
    EXTRA_SOURCES check_reunit_deriv.c
)
