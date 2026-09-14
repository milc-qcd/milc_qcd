# Source groups and targets follow Make_template.
# Legacy machine-dependent I/O uses the current serial/MPI communication
# and ANSI I/O interfaces; old source names are retained below.

add_interface_library(my_objects
    fuz_prop.c
    fuz_source.c
    gauge_info.c
    li_li_mesons.c
    setup.c
    smearing.c
    stat_li_mesons.c
    su3vecsrc_copy.c
    w_loop1.c
    w_loop2_red1.c
    wl_1l_1corr_new.c
    wl_1l_1corr_offaxn_red1_new.c
    wl_1l_2corr_new.c
    wl_1l_2corr_offaxn_red1_new.c
    wl_2l_1corr_new.c
    wl_2l_1corr_offaxn_red1_new.c
    wl_2l_2corr_new.c
    wl_2l_2corr_offaxn_red1_new.c
    # TODO: Restore the source referenced by Make_template.
    layout_hyper.c
)

add_interface_library(g_objects
    ../generic/ape_smear.c
    ../generic/ax_gauge.c
    ../generic/check_unitarity.c
    ../generic_ks/d_congrad5.c
    ../generic_ks/dslash.c
    ../generic/d_plaq4.c
    ../generic/io_helpers.c
    ../generic/io_lat4.c
    ../generic/make_lattice.c
    ../generic/nersc_cksum.c
    ../generic/rand_gauge.c
    ../generic/ranstuff.c
    ../generic/reunitarize2.c
    ../generic/project_su3_hit.c
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

add_milc_executable(su3_str_brk_rgt control5.c
    DEFINES RAN_GAUGE
)

add_milc_executable(su3_str_brk control5.c
)

