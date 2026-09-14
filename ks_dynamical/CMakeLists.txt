# Source groups and targets follow Make_template.
# Legacy machine-dependent I/O uses the current serial/MPI communication
# and ANSI I/O interfaces; old source names are retained below.

add_interface_library(my_objects
    f_measure.c
    gauge_info.c
    setup.c
    update.c
    update_h.c
    update_u.c
)

add_interface_library(g_objects
    ../generic/ranmom.c
    ../generic/layout_squares.c
    ../generic_ks/d_congrad5.c
    ../generic/check_unitarity.c
    ../generic/d_plaq4.c
    ../generic/io_helpers.c
    ../generic/io_lat4.c
    ../generic/make_lattice.c
    ../generic/ploop3.c
    ../generic/ranstuff.c
    ../generic/reunitarize2.c
    ../generic_ks/grsource.c
    ../generic_ks/rephase.c
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
    EXTRA_SOURCES
    analyze.c
)

add_milc_executable(su3_phi control.c
    DEFINES PHI_ALGORITHM
    EXTRA_SOURCES
    analyze.c
)

add_milc_executable(su3_hmc control.c
    DEFINES PHI_ALGORITHM HMC_ALGORITHM
    EXTRA_SOURCES
    d_action.c
    analyze.c
)

add_milc_executable(su3_rmd_spectrum control.c
    DEFINES SPECTRUM
    EXTRA_SOURCES
    ../generic/gaugefix2.c
    ../generic_ks/spectrum.c
    analyze_spectrum.c
)

add_milc_executable(su3_rmd_nlspec control.c
    DEFINES NLSPEC
    EXTRA_SOURCES
    ../generic/gaugefix2.c
    ../generic_ks/spectrum.c
    mom_spec5.c
    ../generic_ks/nl_spectrum.c
    analyze_nlspec.c
)

add_milc_executable(su3_rmd_screen control.c
    DEFINES SPECTRUM
    EXTRA_SOURCES
    ../generic/gaugefix2.c
    spectrum_s2.c
    analyze_spectrum.c
)

add_milc_executable(su3_rmd_screenfunc control.c
    DEFINES SCREENFUNC
    EXTRA_SOURCES
    ../generic/gaugefix2.c
    spectrum_s2.c
    wavefunc_z.c
    analyze_screenfunc.c
)

