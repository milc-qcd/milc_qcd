# Source groups and targets follow Make_template.
# Legacy machine-dependent I/O uses the current serial/MPI communication
# and ANSI I/O interfaces; old source names are retained below.

add_interface_library(my_objects
    bicgstab.c
    build_params.c
    delta0.c
    f_mu_nu1.c
    fp_info.c
    gauge_info.c
    make_clov1.c
    monte_block_ape_b_many.c
    mult_ldu1.c
    path.c
    setup_links.c
    setup_offset.c
    setup_p_cl.c
    # TODO: Restore the source referenced by Make_template.
    layout_hyper.c
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
    ../generic/reunitarize2.c
    ../generic/restrict_fourier.c
    # TODO: Restore the source referenced by Make_template.
    ../generic_wilson/w_source.c
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

add_milc_executable(su3_fp_cl_bi control_p.c
    DEFINES SPECTRUM CVC BI CGTIME
    EXTRA_SOURCES
    cvc_fv.c
    hax.c
    cvc_pauli.c
    w_meson_p.c
    w_baryon_p.c
    rotate_field.c
)

add_milc_executable(su3_fp_cl_bi_h control_p.c
    DEFINES SPECTRUM BI CGTIME
    EXTRA_SOURCES
    w_meson_p.c
    w_baryon_p.c
    rotate_field.c
)

add_milc_executable(su3_fp_cl_bi_n control_p_nothing.c
    DEFINES SPECTRUM BI CGTIME
)

add_milc_executable(su3_fp_cl_eig control_eig.c
    DEFINES RANDOM EIG BI CGTIME
    EXTRA_SOURCES
    ../generic_wilson/boundary_flip.c
    ../generic/ranstuff.c
    f_measure2.c
)

add_milc_executable(su3_fp_cg_eig control_eig.c
    DEFINES RANDOM EIG CGTIME
    EXTRA_SOURCES
    ../generic_wilson/boundary_flip.c
    ../generic/ranstuff.c
    f_measure2.c
    congrad_t.c
)

