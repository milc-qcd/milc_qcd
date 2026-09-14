# Source groups and targets follow Make_template.
# Legacy machine-dependent I/O uses the current serial/MPI communication
# and ANSI I/O interfaces; old source names are retained below.

add_interface_library(my_objects
    setup.c
    gauge_info.c
    # TODO: Restore the source referenced by Make_template.
    layout_hyper.c
)

add_interface_library(g_objects
    ../generic/check_unitarity.c
    ../generic/io_helpers.c
    ../generic/io_lat4.c
    ../generic/make_lattice.c
    ../generic/nersc_cksum.c
    ../generic/plaquette4.c
    ../generic/d_plaq2.c
    ../generic/ranstuff.c
    ../generic/reunitarize2.c
    ../generic_schroed/coupling.c
    ../generic_schroed/make_schroed_lattice.c
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

add_milc_executable(su3_schr_rmd control.c
    DEFINES RMD_ALGORITHM
    EXTRA_SOURCES
    ../generic_pg/update.c
    ../generic_pg/update_u.c
    ../generic_pg/update_h.c
    ../generic/ranmom.c
)

add_milc_executable(su3_schr_hmc control.c
    DEFINES HMC_ALGORITHM
    EXTRA_SOURCES
    ../generic_pg/update.c
    ../generic_pg/update_u.c
    ../generic_pg/update_h.c
    ../generic/ranmom.c
    ../generic_pg/d_action.c
)

add_milc_executable(su3_schr_ora control.c
    DEFINES ORA_ALGORITHM
    EXTRA_SOURCES
    ../generic_pg/update_ora.c
    ../generic_pg/dsdu_qhb.c
    ../generic_pg/relax.c
    ../generic_pg/monte.c
)

