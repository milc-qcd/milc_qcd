# Source groups and targets follow Make_template.

add_interface_library(my_objects
    setup.c
    gauge_info.c
)

add_interface_library(g_objects
    ../generic/check_unitarity.c
    ../generic/d_plaq4.c
    ../generic/field_utilities.c
    ../generic/gauge_action_imp.c
    ../generic/gauge_stuff.c
    ../generic/gauge_force_imp.c
    ../generic/gauge_utilities.c
    ../generic/io_detect.c
    ../generic/io_helpers.c
    ../generic/io_lat4.c
    ../generic/io_lat_utils.c
    ../generic/make_lattice.c
    ../generic/nersc_cksum.c
    ../generic/remap_stdio_from_args.c
    ../generic/path_product.c
    ../generic/plaquette4.c
    ../generic/ploop2.c
    ../generic/ranstuff.c
    ../generic/reunitarize2.c
    ../generic/layout_hyper_sl32.c
)

if(HAVEQIO)
    append_interface_library(g_objects
        ../generic/file_types_milc_usqcd.c
        ../generic/io_scidac.c
        ../generic/io_scidac_types.c
        ../generic_ks/io_scidac_ks.c
        ../generic_wilson/io_scidac_w.c
    )
endif()

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

add_milc_executable(su3_hmc_symzk1 control.c
    GAUGE symanzik_1loop_action
    DEFINES HMC_ALGORITHM
    EXTRA_SOURCES
    update.c
    update_u.c
    update_h.c
    ../generic/ranmom.c
    d_action.c
)

add_milc_executable(su3_ora_symzk1 control.c
    GAUGE symanzik_1loop_action
    DEFINES ORA_ALGORITHM
    EXTRA_SOURCES
    update_ora.c
    relax.c
    monte.c
)

