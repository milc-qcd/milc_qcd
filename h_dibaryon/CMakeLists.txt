# Source groups and targets follow Make_template.
# Legacy machine-dependent I/O uses the current serial/MPI communication
# and ANSI I/O interfaces; old source names are retained below.

add_interface_library(my_objects
    clover_info.c
    diquarkprop.c
    gauge_info.c
    make_channel_wfs.c
    nr_propagator.c
    setup_H_cl.c
    w_hdibaryon.c
    w_nrbaryon.c
    # TODO: Restore the source referenced by Make_template.
    layout_hyper.c
    # TODO: Restore the source referenced by Make_template.
    dslash_lean.c
    # TODO: Restore the source referenced by Make_template.
    wilson_invert_lean.c
    # TODO: Restore the source referenced by Make_template.
    make_clov.c
)

add_interface_library(g_objects
    ../generic/check_unitarity.c
    ../generic/d_plaq4.c
    ../generic/gaugefix2.c
    ../generic/gauge_utilities.c
    ../generic/io_helpers.c
    ../generic_wilson/io_helpers_w.c
    ../generic/io_lat4.c
    ../generic_wilson/io_prop_w.c
    ../generic/make_lattice.c
    ../generic/restrict_fourier.c
    ../generic/reunitarize2.c
    ../generic_wilson/meson_cont.c
    ../generic_wilson/w_meson.c
    # TODO: Restore the source referenced by Make_template.
    ../generic_wilson/w_source.c
    ../generic_clover/f_mu_nu.c
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

add_milc_executable(su3_Hcl_cg control_H_cl.c
    EXTRA_SOURCES
    ../generic_clover/d_cgilu_cl.c
)

add_milc_executable(su3_Hcl_bi control_H_cl.c
    DEFINES BI
    EXTRA_SOURCES
    ../generic_clover/d_bicgilu_cl.c
)

