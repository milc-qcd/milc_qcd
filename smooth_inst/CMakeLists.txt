add_interface_library(my_objects
    char_num.c
    gauge_info.c
    make_loop_table.c
    path.c
    ../generic/project_su3_hit.c
    save_topo.c
    setup.c
)

add_interface_library(g_objects
    ../generic/check_unitarity.c
    ../generic/d_plaq4.c
    ../generic/gauge_utilities.c
    ../generic/io_helpers.c
    ../generic/io_lat4.c
    ../generic/io_lat_utils.c
    ../generic/make_lattice.c
    ../generic/nersc_cksum.c
    ../generic/remap_stdio_from_args.c
    ../generic/reunitarize2.c
)

add_interface_library(ks_g_objects
    ../generic_ks/fermion_links_hyp.c
)

if(HAVEQIO)
    append_interface_library(g_objects
        ../generic/io_scidac.c
        ../generic/io_scidac_types.c
    )
endif()

add_interface_library(layout
    ../generic/layout_hyper_prime.c
)

add_interface_library(add_objects
    OBJECTS
    machine_dep_io
    commpkg
    ks_g_objects
)

add_interface_library(objects
    OBJECTS
    my_objects
    g_objects
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

add_milc_executable(su3_hyp control.c
    DEFINES HYP AB_USE_GENERIC_HYP
    EXTRA_SOURCES hyp_block.c ../generic/gaugefix2.c instanton_density.c
)

add_milc_executable(su3_hyp2 control.c
    DEFINES HYP
    EXTRA_SOURCES hyp_block.c ../generic/gaugefix2.c instanton_density2.c
)

add_milc_executable(su3_ape control.c
    EXTRA_SOURCES smooth.c ../generic/ape_smear.c ../generic/gaugefix2.c instanton_density.c
)

add_milc_executable(su3_fat control.c
    EXTRA_SOURCES fat.c ../generic/ape_smear.c ../generic/gaugefix2.c instanton_density.c
)

add_milc_executable(su3_stout control.c
    EXTRA_SOURCES smooth_stout.c ../generic/ape_smear.c ../generic/stout_smear.c ../generic/gaugefix2.c instanton_density.c
)
