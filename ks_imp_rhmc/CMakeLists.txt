set(my_sources
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

set(layout_sources
    ../generic/layout_hyper_prime.c
)

set(add_sources
    ../generic/ploop3.c
    ../generic_ks/show_generic_ks_md_opts.c
)

set(sources
    ${my_sources}
    ${layout_sources}
    ${add_sources}
)

add_library(libraries INTERFACE)
target_link_libraries(libraries INTERFACE
    gauge_iface ${gauge_force} standard_iface
    qcd scidac add machine_dep_io comm m
)

add_executable(su3_leapfrog ${sources} control.c)
target_link_action(su3_leapfrog symanzik_1loop_action asqtad_action)
target_compile_definitions(su3_leapfrog PRIVATE FN INT_ALG=INT_LEAPFROG)
target_link_libraries(su3_leapfrog PRIVATE fn_iface ${asq_force} libraries)

add_executable(su3_omelyan_rhmc ${sources} control.c)
target_link_action(su3_omelyan_rhmc symanzik_1loop_action asqtad_action)
target_compile_definitions(su3_omelyan_rhmc PRIVATE FN HMC INT_ALG=INT_OMELYAN)
target_link_libraries(su3_omelyan_rhmc PRIVATE fn_iface ${asq_force} libraries)

set(hisq_options
    FN
    INT_ALG=INT_3G1F
    HISQ_REUNIT_ALLOW_SVD
    HISQ_REUNIT_SVD_REL_ERROR=1e-8
    DHISQ_REUNIT_SVD_ABS_ERROR=1e-8
    DHISQ_FF_MULTI_WRAPPER
    DHISQ_FORCE_FILTER=5.0e-5
)

add_executable(su3_rhmd_hisq ${sources} control.c)
target_link_action(su3_rhmd_hisq symanzik_1loop_hisq_action hisq/hisq_u3_action)
target_compile_definitions(su3_rhmd_hisq PRIVATE ${hisq_options})
target_link_libraries(su3_rhmd_hisq PRIVATE hisq_iface ${hisq_force} libraries)
