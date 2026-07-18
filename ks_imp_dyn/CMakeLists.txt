add_interface_library(my_objects
    gauge_info.c
    ../generic_ks/reunitarize_ks.c
    setup.c
    update_h.c
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

add_interface_library(action_objects
    d_action.c
)

set(ks_multicg_off
    -UKS_MULTICG
)

set(ks_multicg_multiff_off
    -UKS_MULTICG
    -UKS_MULTIFF
)

set(asqtad_base
    FN
)

add_milc_executable(su3_rmd control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES FN
    EXTRA_OBJECTS fn_objects asq_force gauge_force
    EXTRA_SOURCES update.c
)
target_compile_options(su3_rmd PRIVATE ${ks_multicg_off})

add_milc_executable(su3_rmd1 control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES FN ONEMASS
    EXTRA_OBJECTS fn_objects asq_force gauge_force
    EXTRA_SOURCES update_onemass.c
)
target_compile_options(su3_rmd1 PRIVATE ${ks_multicg_off})

add_milc_executable(su3_rmd_eos control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES FN DM_DU0 TR_MM_INV
    EXTRA_OBJECTS fn_objects asq_force gauge_force
    EXTRA_SOURCES update.c
)
target_compile_options(su3_rmd_eos PRIVATE ${ks_multicg_off})

add_milc_executable(su3_rmd_eos1 control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES FN DM_DU0 ONEMASS
    EXTRA_OBJECTS fn_objects asq_force gauge_force
    EXTRA_SOURCES update_onemass.c
)
target_compile_options(su3_rmd_eos1 PRIVATE ${ks_multicg_off})

add_milc_executable(su3_rmd_susc control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES FN CHEM_POT TR_MM_INV
    EXTRA_OBJECTS fn_objects asq_force gauge_force
    EXTRA_SOURCES update.c
)
target_compile_options(su3_rmd_susc PRIVATE ${ks_multicg_off})

add_milc_executable(su3_rmd_mu_eos control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES FN DM_DU0 D_CHEM_POT TR_MM_INV
    EXTRA_OBJECTS fn_objects asq_force gauge_force
    EXTRA_SOURCES update.c ../generic_ks/mu_fast.c
)
target_compile_options(su3_rmd_mu_eos PRIVATE ${ks_multicg_off})

add_milc_executable(su3_rmd_dblstore_gath13_fewsums control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES FN DBLSTORE_FN D_FN_GATHER13 FEWSUMS
    EXTRA_OBJECTS fn_objects asq_force gauge_force
    EXTRA_SOURCES update.c
)
target_compile_options(su3_rmd_dblstore_gath13_fewsums PRIVATE ${ks_multicg_off})

add_milc_executable(su3_hmc control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES FN PHI_ALGORITHM HMC_ALGORITHM
    EXTRA_OBJECTS fn_objects asq_force gauge_force
    EXTRA_SOURCES d_action.c update.c
)
target_compile_options(su3_hmc PRIVATE ${ks_multicg_off})

add_milc_executable(su3_hmc1 control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES FN ONEMASS PHI_ALGORITHM HMC_ALGORITHM
    EXTRA_OBJECTS fn_objects asq_force gauge_force action_objects
    EXTRA_SOURCES update_onemass.c
)
target_compile_options(su3_hmc1 PRIVATE ${ks_multicg_off})

add_milc_executable(su3_rmd_naive control.c
    GAUGE plaquette_action
    QUARK onelink_action
    DEFINES FN ONEMASS
    EXTRA_OBJECTS fn_objects asq_force gauge_force
    EXTRA_SOURCES update_onemass.c
)
target_compile_options(su3_rmd_naive PRIVATE ${ks_multicg_multiff_off} -UD_FN_GATHER13)

add_milc_executable(su3_hmc_eo_symzk1_p4 control.c
    GAUGE symanzik_1loop_action
    QUARK p4_action
    DEFINES EO PHI_ALGORITHM HMC_ALGORITHM
    EXTRA_OBJECTS eo_objects eo_force gauge_force_milc
    EXTRA_SOURCES ../generic_ks/d_congrad5_two_src.c d_action.c update.c
)
target_compile_options(su3_hmc_eo_symzk1_p4 PRIVATE ${ks_multicg_multiff_off} -UHAVE_QOP -UHAVE_QDP)

add_milc_executable(su3_rmd_eo_symzk1_fat7tad control.c
    GAUGE symanzik_1loop_action
    QUARK fat7tad_action
    DEFINES EO
    EXTRA_OBJECTS eo_objects gauge_force
    EXTRA_SOURCES update.c
)
target_compile_options(su3_rmd_eo_symzk1_fat7tad PRIVATE ${ks_multicg_off})

add_milc_executable(su3_rmd_plaq_fat7naik control.c
    GAUGE plaquette_action
    QUARK fat7naik_action
    DEFINES FN
    EXTRA_OBJECTS fn_objects asq_force gauge_force
    EXTRA_SOURCES update.c
)
target_compile_options(su3_rmd_plaq_fat7naik PRIVATE ${ks_multicg_multiff_off})

add_milc_executable(su3_rmd_plaq_fat7naik1 control.c
    GAUGE plaquette_action
    QUARK fat7naik_action
    DEFINES FN ONEMASS
    EXTRA_OBJECTS fn_objects asq_force gauge_force
    EXTRA_SOURCES update_onemass.c
)
target_compile_options(su3_rmd_plaq_fat7naik1 PRIVATE ${ks_multicg_multiff_off})

add_milc_executable(su3_rmd_symzk1_naik control.c
    GAUGE symanzik_1loop_action
    QUARK naik_action
    DEFINES FN
    EXTRA_OBJECTS fn_objects asq_force gauge_force
    EXTRA_SOURCES update.c
)
target_compile_options(su3_rmd_symzk1_naik PRIVATE ${ks_multicg_off})

add_milc_executable(su3_rmd_eo_symzk1_p4 control.c
    GAUGE symanzik_1loop_action
    QUARK p4_action
    DEFINES EO
    EXTRA_OBJECTS eo_objects eo_force gauge_force
    EXTRA_SOURCES update.c
)
target_compile_options(su3_rmd_eo_symzk1_p4 PRIVATE ${ks_multicg_multiff_off} -UHAVE_QOP -UHAVE_QDP)

set(hisq_options
    FN
    HISQ_REUNIT_ALLOW_SVD
    HISQ_REUNIT_SVD_REL_ERROR=1e-8
    HISQ_REUNIT_SVD_ABS_ERROR=1e-8
    HISQ_FF_MULTI_WRAPPER
    HISQ_FORCE_FILTER=5.0e-5
)

set(hisq_tree_options
    FN
    HISQ_REUNIT_ALLOW_SVD
    HISQ_REUNIT_SVD_REL_ERROR=1e-8
    HISQ_REUNIT_SVD_ABS_ERROR=1e-8
    HISQ_FF_MULTI_WRAPPER
    HISQ_FORCE_FILTER=5.0e-5
)

add_milc_executable(su3_rmd_hisq control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options}
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
    EXTRA_SOURCES update.c
)
target_compile_options(su3_rmd_hisq PRIVATE ${ks_multicg_off})

add_milc_executable(su3_rmd_hisq_debug control.c
    GAUGE symanzik_1loop_action
    QUARK hisq/asqtad_action
    DEFINES ${hisq_options}
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
    EXTRA_SOURCES update.c
)
target_compile_options(su3_rmd_hisq_debug PRIVATE ${ks_multicg_off})

add_milc_executable(su3_rmd1_hisq control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} ONEMASS
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
    EXTRA_SOURCES update_onemass.c
)
target_compile_options(su3_rmd1_hisq PRIVATE ${ks_multicg_off})

add_milc_executable(su3_hmc_hisq control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} HMC_ALGORITHM
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
    EXTRA_SOURCES d_action.c update.c
)
target_compile_options(su3_hmc_hisq PRIVATE ${ks_multicg_off})

add_milc_executable(su3_hmc1_hisq control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} ONEMASS PHI_ALGORITHM HMC_ALGORITHM
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
    EXTRA_SOURCES d_action.c update_onemass.c
)
target_compile_options(su3_hmc1_hisq PRIVATE ${ks_multicg_off})

add_milc_executable(su3_hmc1_tree_hisq control.c
    GAUGE symanzik_tree_action
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_tree_options}
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
    EXTRA_SOURCES d_action.c update_onemass.c
)
target_compile_options(su3_hmc1_tree_hisq PRIVATE ${ks_multicg_off})
