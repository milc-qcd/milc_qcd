add_interface_library(my_objects
    gauge_action.c
    loop_table.c
    mult_sigma_mu_nu.c
    path.c
    setup.c
    udadu_mat_mu_nu.c
    udadu_mu_nu.c
    update.c
    update_u.c
    update_h_cl.c
    grsource.c
    f_measure_cl.c
    gauge_info.c
)

add_interface_library(g_objects
    ../generic/check_unitarity.c
    ../generic/d_plaq4.c
    ../generic/field_utilities.c
    ../generic/io_helpers.c
    ../generic/io_lat4.c
    ../generic/io_lat_utils.c
    ../generic/make_lattice.c
    ../generic/nersc_cksum.c
    ../generic/ploop_dist.c
    ../generic/ranmom.c
    ../generic/ranstuff.c
    ../generic/remap_stdio_from_args.c
    ../generic/report_invert_status.c
    ../generic/reunitarize2.c
)
if(HAVEQIO)
    append_interface_library(g_objects
        ../generic/io_scidac.c
        ../generic/io_scidac_types.c
    )
endif()

add_interface_library(g_wilson_objects
    ../generic_wilson/boundary_flip.c
    ../generic_wilson/wilson_invert.c
    ../generic_wilson/dslash_w3.c
)

add_interface_library(g_clover_objects
    ../generic_clover/cl_solver_utilities.c
    ../generic_clover/d_congrad2_cl.c
    ../generic_clover/d_cgilu_cl.c
    ../generic_clover/f_mu_nu.c
    ../generic_clover/make_clov2.c
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
    g_objects
    g_wilson_objects
    g_clover_objects
    layout
    add_objects
)

add_interface_library(libraries
    OBJECTS
    qcd
    scidac
    add
    m
)

add_milc_executable(su3_rmd control.c
    DEFINES LU
)

add_milc_executable(su3_rmd_bi control.c
    DEFINES LU BI
    EXTRA_SOURCES ../generic_clover/d_bicgilu_cl.c
)

add_milc_executable(su3_phi control.c
    DEFINES LU PHI_ALGORITHM BI
)

add_milc_executable(su3_hmc control.c
    DEFINES LU PHI_ALGORITHM HMC_ALGORITHM
    EXTRA_SOURCES d_action.c
)

add_milc_executable(su3_hmc_bi control.c
    DEFINES LU PHI_ALGORITHM HMC_ALGORITHM BI
    EXTRA_SOURCES d_action.c ../generic_clover/d_bicgilu_cl.c
)

add_milc_executable(su3_hmc_spectrum control.c
    DEFINES LU PHI_ALGORITHM HMC_ALGORITHM SPECTRUM
    EXTRA_SOURCES
    d_action.c
    ../generic/gaugefix2.c
    t_props_cl.c
    w_spectrum_cl.c
)

add_milc_executable(su3_hmc_bi_spectrum control.c
    DEFINES LU PHI_ALGORITHM HMC_ALGORITHM BI SPECTRUM
    EXTRA_SOURCES
    d_action.c
    ../generic/gaugefix2.c
    t_props_cl.c
    w_spectrum_cl.c
    ../generic_clover/d_bicgilu_cl.c
)

add_milc_executable(su3_rmd_spectrum control.c
    DEFINES LU SPECTRUM
    EXTRA_SOURCES
    ../generic/gaugefix2.c
    t_props_cl.c
    ../wilson_dynamical/w_spectrum.c
)

add_milc_executable(su3_hmc_screen control.c
    DEFINES LU PHI_ALGORITHM HMC_ALGORITHM SPECTRUM SCREEN
    EXTRA_SOURCES
    d_action.c
    ../generic/gaugefix2.c
    s_props_cl.c
)

add_milc_executable(su3_hmc_bi_screen control.c
    DEFINES LU PHI_ALGORITHM HMC_ALGORITHM BI SPECTRUM SCREEN
    EXTRA_SOURCES
    d_action.c
    ../generic/gaugefix2.c
    s_props_cl.c
    ../generic_clover/d_bicgilu_cl.c
)
