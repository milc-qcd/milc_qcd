# Source groups and targets follow Make_template.
# Legacy machine-dependent I/O uses the current serial/MPI communication
# and ANSI I/O interfaces; old source names are retained below.

add_interface_library(my_objects
    # TODO: Restore the source referenced by Make_template.
    d_congrad2.c
    setup.c
    update.c
    update_u.c
    update_h.c
    grsource_w.c
    f_measure3.c
    gauge_info.c
    # TODO: Restore the source referenced by Make_template.
    layout_hyper.c
    # TODO: Restore the source referenced by Make_template.
    dslash_w.c
)

add_interface_library(g_objects
    ../generic/ranmom.c
    ../generic/d_plaq4.c
    ../generic/field_utilities.c
    ../generic/io_lat4.c
    ../generic/io_helpers.c
    ../generic/check_unitarity.c
    ../generic/make_lattice.c
    ../generic/nersc_cksum.c
    ../generic/reunitarize2.c
    ../generic/ploop_dist.c
    ../generic/ranstuff.c
    ../generic_wilson/boundary_flip.c
)

add_interface_library(objects
    OBJECTS
    my_objects
    g_objects
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

add_milc_executable(su3_rmd control.c
    DEFINES LU
)

add_milc_executable(su3_phi control.c
    DEFINES LU PHI_ALGORITHM
)

add_milc_executable(su3_hmc control.c
    DEFINES LU PHI_ALGORITHM HMC_ALGORITHM
    EXTRA_SOURCES
    d_action.c
)

add_milc_executable(su3_hmc_spectrum control.c
    DEFINES LU PHI_ALGORITHM HMC_ALGORITHM SPECTRUM
    EXTRA_SOURCES
    d_action.c
    ../generic/gaugefix2.c
    t_props.c
    w_spectrum.c
)

add_milc_executable(su3_hmc_screen control.c
    DEFINES LU PHI_ALGORITHM HMC_ALGORITHM SPECTRUM SCREEN
    EXTRA_SOURCES
    d_action.c
    ../generic/gaugefix2.c
    s_props.c
    w_spectrum_z2.c
)

add_milc_executable(su3_rmd_spectrum control.c
    DEFINES LU SPECTRUM
    EXTRA_SOURCES
    ../generic/gaugefix2.c
    t_props.c
    w_spectrum.c
)

add_milc_executable(su3_rmd_screen control.c
    DEFINES LU SPECTRUM SCREEN
    EXTRA_SOURCES
    d_action.c
    ../generic/gaugefix2.c
    s_props.c
    w_spectrum_z2.c
)

