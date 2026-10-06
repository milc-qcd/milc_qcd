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

add_milc_executable(su3_omelyan_rhmc control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES FN HMC INT_ALG=INT_OMELYAN
    EXTRA_OBJECTS fn_objects asq_force gauge_force
)

set(asq_options_nonnf
    FN
)

add_interface_library(debug_objects
    debug.c
)

add_milc_executable(su3_rhmd control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES ${asq_options_nonnf}
    EXTRA_OBJECTS fn_objects asq_force gauge_force
)

add_milc_executable(su3_rhmc control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES ${asq_options_nonnf} HMC
    EXTRA_OBJECTS fn_objects asq_force gauge_force
)

add_milc_executable(su3_rhmd_eos control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES ${asq_options_nonnf} DM_DU0
    EXTRA_OBJECTS fn_objects asq_force gauge_force
)

add_milc_executable(su3_rhmc_eos control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES ${asq_options_nonnf} HMC DM_DU0
    EXTRA_OBJECTS fn_objects asq_force gauge_force
)

add_milc_executable(su3_rhmc_susc control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES ${asq_options_nonnf} HMC CHEM_POT
    EXTRA_OBJECTS fn_objects asq_force gauge_force
)

add_milc_executable(su3_rhmc_susc_eos control.c
    GAUGE symanzik_1loop_action
    QUARK asqtad_action
    DEFINES ${asq_options_nonnf} HMC CHEM_POT DM_DU0
    EXTRA_OBJECTS fn_objects asq_force gauge_force
)

set(asq_options_nf
    FN
)

add_milc_executable(su3_rhmc_Nf control.c
    GAUGE symanzik_1loop_nf_action
    QUARK asqtad_action
    DEFINES ${asq_options_nf} HMC
    EXTRA_OBJECTS fn_objects asq_force gauge_force
)

set(hisq_options
    FN
    INT_ALG=INT_3G1F
    HISQ_REUNIT_ALLOW_SVD
    HISQ_REUNIT_SVD_REL_ERROR=1e-8
    HISQ_REUNIT_SVD_ABS_ERROR=1e-8
    HISQ_FF_MULTI_WRAPPER
    HISQ_FORCE_FILTER=5.0e-5
)

add_milc_executable(su3_rhmd_hisq control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options}
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)

add_milc_executable(su3_rhmc_hisq control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} HMC
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)

add_milc_executable(su3_rhmd_hisq_su3 control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_su3_action
    DEFINES FN INT_ALG=INT_3G1F
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)

add_milc_executable(su3_rhmc_hisq_5G1F control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_u3_action
    DEFINES FN HMC INT_ALG=INT_5G1F HISQ_REUNIT_ALLOW_SVD HISQ_REUNIT_SVD_REL_ERROR=1e-8 HISQ_REUNIT_SVD_ABS_ERROR=1e-8 HISQ_FF_MULTI_WRAPPER HISQ_FORCE_FILTER=5.0e-5
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)

add_milc_executable(su3_rhmc_hisq_6G1F control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_u3_action
    DEFINES FN HMC INT_ALG=INT_6G1F HISQ_REUNIT_ALLOW_SVD HISQ_REUNIT_SVD_REL_ERROR=1e-8 HISQ_REUNIT_SVD_ABS_ERROR=1e-8 HISQ_FF_MULTI_WRAPPER HISQ_FORCE_FILTER=5.0e-5
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)

add_milc_executable(su3_rhmc_hisq_plaq control.c
    GAUGE symanzik_1loop_nf0_action
    QUARK hisq/hisq_u3_action
    DEFINES FN HMC INT_ALG=INT_3G1F HISQ_REUNIT_ALLOW_SVD HISQ_REUNIT_SVD_ONLY
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)

add_milc_executable(su3_rhmc_hisq_su3 control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_su3_action
    DEFINES FN HMC INT_ALG=INT_3G1F
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)

add_milc_executable(su3_rhmc_hisq_eos control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} HMC DM_DEPS0
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)

add_milc_executable(su3_rhmd_hisq_eos control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} DM_DEPS
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)

add_milc_executable(su3_rhmc_hisq_susc control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} HMC CHEM_POT
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)

add_milc_executable(su3_rhmc_hisq_susc_eos control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} HMC CHEM_POT DM_DEPS
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)

add_milc_executable(su3_rhmc_asq_debug control.c
    GAUGE symanzik_1loop_action
    QUARK asq_action
    DEFINES FN HMC INT_ALG=INT_3G1F MILC_GLOBAL_DEBUG
    EXTRA_OBJECTS fn_objects asq_force gauge_force debug_objects
)

add_milc_executable(su3_rhmc_hisq_debug control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_u3_action
    DEFINES FN HMC INT_ALG=INT_3G1F MILC_GLOBAL_DEBUG HISQ_REUNITARIZATION_DEBUG HISQ_FF_MULTI_WRAPPERRR HISQ_FF_DEBUGGG
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force debug_objects
)

add_milc_executable(su3_rhmc_hisq_su3_debug control.c
    GAUGE symanzik_1loop_hisq_action
    QUARK hisq/hisq_su3_action
    DEFINES FN HMC INT_ALG=INT_3G1F MILC_GLOBAL_DEBUG HISQ_REUNITARIZATION_DEBUG HISQ_FF_MULTI_WRAPPER HISQ_FF_DEBUG
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force debug_objects
)

add_milc_executable(su3_rhmc_hisq_wrap_asq control.c
    GAUGE symanzik_1loop_action
    QUARK hisq/1Xasqtad_action
    DEFINES FN HMC INT_ALG=INT_3G1F MILC_GLOBAL_DEBUG HISQ_FF_MULTI_WRAPPER
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force debug_objects
)

add_milc_executable(su3_rhmc_hisq_tune control.c
    GAUGE symanzik_1loop_nf_action
    QUARK hisq/hisq_u3_action
    DEFINES FN HMC MILC_GLOBAL_DEBUG MEASURE_AND_TUNE_HISQ
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force debug_objects
)

add_milc_executable(test_su3_mat_op_inverse control_test_su3_mat_op.c
    GAUGE symanzik_1loop_nf_action
    QUARK hisq/hisq_action
    DEFINES FN TEST_SU3_INVERSION
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)

add_milc_executable(test_su3_mat_op_unit control_test_su3_mat_op.c
    GAUGE symanzik_1loop_nf_action
    QUARK hisq/hisq_u3_action
    DEFINES FN HISQ_REUNIT_ALLOW_SVD HISQ_REUNIT_SVD_ONLY TEST_SU3_UNITARIZATION
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)

add_milc_executable(test_su3_mat_op_unit_der control_test_su3_mat_op.c
    GAUGE symanzik_1loop_nf_action
    QUARK hisq/hisq_u3_action
    DEFINES FN TEST_SU3_UNITARIZATION_DERIVATIVE
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)

add_milc_executable(test_su3_paranoid control_test_su3_mat_op.c
    GAUGE symanzik_1loop_nf_action
    QUARK hisq/hisq_u3_action
    DEFINES FN HISQ_REUNIT_ALLOW_SVD TEST_SU3_UNITARIZATION
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)

add_milc_executable(test_su3_paranoid_der control_test_su3_mat_op.c
    GAUGE symanzik_1loop_nf_action
    QUARK hisq/hisq_u3_action
    DEFINES FN HISQ_REUNIT_ALLOW_SVD TEST_SU3_UNITARIZATION_DERIVATIVE
    EXTRA_OBJECTS hisq_objects hisq_force gauge_force
)
