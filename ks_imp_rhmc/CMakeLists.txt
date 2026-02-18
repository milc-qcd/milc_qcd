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

add_library(library INTERFACE)
target_link_libraries(library INTERFACE qcd scidac add machine_dep_io comm)

set(hisq_options
    FN INT_ALG=INT_3G1F
    HISQ_REUNIT_ALLOW_SVD
    HISQ_REUNIT_SVD_REL_ERROR=1e-8
    DHISQ_REUNIT_SVD_ABS_ERROR=1e-8
    DHISQ_FF_MULTI_WRAPPER
    DHISQ_FORCE_FILTER=5.0e-5
)

add_executable(su3_rhmd_hisq ${sources} control.c)
target_link_action(su3_rhmd_hisq symanzik_1loop_hisq_action hisq/hisq_u3_action)
target_link_libraries(su3_rhmd_hisq PRIVATE hisq_iface ${hisq_force} gauge_iface ${gauge_force} standard_iface library)
target_compile_definitions(su3_rhmd_hisq PRIVATE ${hisq_options})
target_include_directories(su3_rhmd_hisq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
