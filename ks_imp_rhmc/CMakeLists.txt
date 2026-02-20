add_interface_library(my_objects
    d_action_rhmc.c
    eo_fermion_force_rhmc.c
    ../generic_ks/gauge_action_imp_ks.c
    gauge_info.c
    ../generic_ks/gauge_force_imp_ks.c
    ../generic_ks/gauge_measure_ks.c
    ../generic_ks/grsource_rhmc.c
    ../generic_ks/ks_ratinv.c
    ../generic_ks/load_rhmc_params.c
    ../generic_ks/reunitarize_ks.c
    ../generic_ks/ploop3_ks.c
    setup.c
    update_h_rhmc.c
    update_rhmc.c
    update_u.c
)

add_interface_library(layout
    ../generic/layout_hyper_prime.c
)

add_interface_library(add_objects
    ../generic/ploop3.c
    ../generic_ks/show_generic_ks_md_opts.c
    OBJECTS
    machine_dep_io
    commpkg
)

add_interface_library(objects
    OBJECTS
    my_objects
    gauge_objects
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

add_milc_executable(su3_leapfrog control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES FN INT_ALG=INT_LEAPFROG
    EXTRA_OBJECTS fn_objects asq_force gauge_force
)

add_milc_executable(su3_omelyan control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES FN HMC INT_ALG=INT_OMELYAN
    EXTRA_OBJECTS fn_objects asq_force gauge_force
)

set(hisq_options
    FN
    INT_ALG=INT_3G1F
    HISQ_REUNIT_ALLOW_SVD
    HISQ_REUNIT_SVD_REL_ERROR=1e-8
    DHISQ_REUNIT_SVD_ABS_ERROR=1e-8
    DHISQ_FF_MULTI_WRAPPER
    DHISQ_FORCE_FILTER=5.0e-5
)

add_milc_executable(su3_rhmd_hisq control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options}
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)
