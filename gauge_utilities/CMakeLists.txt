add_interface_library(my_objects
    gauge_info.c
    setup.c
)

add_interface_library(g_objects
    ../generic/check_unitarity.c
    ../generic/d_plaq4.c
    ../generic/field_translation.c
    ../generic/field_utilities.c
    ../generic/gauge_utilities.c
    ../generic/gaugefix2.c
    ../generic/gauge_stuff.c
    ../generic/io_detect.c
    ../generic/io_helpers.c
    ../generic/io_lat4.c
    ../generic/io_lat_utils.c
    ../generic/make_lattice.c
    ../generic/momentum_twist.c
    ../generic/nersc_cksum.c
    ../generic/path_product.c
    ../generic/ploop3.c
    ../generic/remap_stdio_from_args.c
    ../generic/reunitarize2.c
    ../generic/ranstuff.c
)
if(HAVEQIO)
    append_interface_library(g_objects
        ../generic/file_types_milc_usqcd.c
        ../generic/io_scidac.c
        ../generic/io_scidac_types.c
    )
endif()

add_interface_library(layout
    ../generic/layout_hyper_prime.c
)

add_interface_library(g_ks_objects
)

add_interface_library(g_w_objects
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
    g_ks_objects
    g_w_objects
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

add_milc_executable(gauge_modify_su3 control.c
    GAUGE symanzik_1loop_action
)
