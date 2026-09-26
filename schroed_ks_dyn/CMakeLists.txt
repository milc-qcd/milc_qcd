# Source groups and targets follow Make_template.

add_interface_library(my_objects
    f_measure.c
    gauge_info.c
    setup.c
    update_u.c
    update_h.c
)

add_interface_library(g_objects
    ../generic/plaquette4.c
    ../generic/ranmom.c
    ../generic_schroed/coupling.c
    ../generic_schroed/make_schroed_lattice.c
    ../generic/layout_hyper_prime.c
)

add_interface_library(objects
    OBJECTS
    my_objects
    g_objects
    gauge_objects
    standard_objects
    fn_objects
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

add_milc_executable(su3_schr_rmd control.c
    QUARK onelink_action
    DEFINES FN
    EXTRA_SOURCES
    update.c
)
target_compile_options(su3_schr_rmd PRIVATE -UD_FN_GATHER13)

add_milc_executable(su3_schr_phi control.c
    QUARK onelink_action
    DEFINES FN PHI_ALGORITHM
    EXTRA_SOURCES
    update.c
)
target_compile_options(su3_schr_phi PRIVATE -UD_FN_GATHER13)

add_milc_executable(su3_schr_hmc control.c
    QUARK onelink_action
    DEFINES FN PHI_ALGORITHM HMC_ALGORITHM
    EXTRA_SOURCES
    d_action.c
    update.c
)
target_compile_options(su3_schr_hmc PRIVATE -UD_FN_GATHER13)

add_milc_executable(su3_schr_sw_phi control.c
    QUARK onelink_action
    DEFINES FN PHI_ALGORITHM SEXT_WEIN
    EXTRA_SOURCES
    update_hgh.c
)
target_compile_options(su3_schr_sw_phi PRIVATE -UD_FN_GATHER13)

add_milc_executable(su3_schr_sw_hmc control.c
    QUARK onelink_action
    DEFINES FN PHI_ALGORITHM HMC_ALGORITHM SEXT_WEIN
    EXTRA_SOURCES
    d_action.c
    update_hgh.c
)
target_compile_options(su3_schr_sw_hmc PRIVATE -UD_FN_GATHER13)

add_milc_executable(su3_schr_rmd_ph control.c
    QUARK onelink_action
    DEFINES FN FERM_PHASES
    EXTRA_SOURCES
    update.c
)
target_compile_options(su3_schr_rmd_ph PRIVATE -UD_FN_GATHER13)

add_milc_executable(su3_schr_phi_ph control.c
    QUARK onelink_action
    DEFINES FN FERM_PHASES PHI_ALGORITHM
    EXTRA_SOURCES
    update.c
)
target_compile_options(su3_schr_phi_ph PRIVATE -UD_FN_GATHER13)

add_milc_executable(su3_schr_hmc_ph control.c
    QUARK onelink_action
    DEFINES FN FERM_PHASES PHI_ALGORITHM HMC_ALGORITHM
    EXTRA_SOURCES
    d_action.c
    update.c
)
target_compile_options(su3_schr_hmc_ph PRIVATE -UD_FN_GATHER13)

add_milc_executable(su3_schr_rw_rmd control.c
    QUARK onelink_action
    DEFINES FN REWEIGH
    EXTRA_SOURCES
    update.c
)
target_compile_options(su3_schr_rw_rmd PRIVATE -UD_FN_GATHER13)

add_milc_executable(su3_schr_rw_phi control.c
    QUARK onelink_action
    DEFINES FN REWEIGH PHI_ALGORITHM
    EXTRA_SOURCES
    update.c
)
target_compile_options(su3_schr_rw_phi PRIVATE -UD_FN_GATHER13)

add_milc_executable(su3_schr_rw_hmc control.c
    QUARK onelink_action
    DEFINES FN REWEIGH PHI_ALGORITHM HMC_ALGORITHM
    EXTRA_SOURCES
    d_action.c
    update.c
)
target_compile_options(su3_schr_rw_hmc PRIVATE -UD_FN_GATHER13)

