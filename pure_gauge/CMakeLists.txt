add_interface_library(my_objects
    setup.c
    gauge_info.c
)

add_interface_library(g_objects
    ../generic/check_unitarity.c
    ../generic/d_plaq4.c
    ../generic/field_utilities.c
    ../generic/io_detect.c
    ../generic/io_helpers.c
    ../generic/io_lat4.c
    ../generic/io_lat_utils.c
    ../generic/make_lattice.c
    ../generic/nersc_cksum.c
    ../generic/ploop2.c
    ../generic/remap_stdio_from_args.c
    ../generic/reunitarize2.c
    ../generic/ranstuff.c
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

add_milc_executable(su3_rmd control.c
    DEFINES RMD_ALGORITHM ONLY_GLUON_FILES
    EXTRA_SOURCES
    ../generic_pg/update.c
    ../generic_pg/update_u.c
    ../generic_pg/update_h.c
    ../generic/ranmom.c
)

add_milc_executable(su3_rmd_a control.c
    DEFINES RMD_ALGORITHM ONLY_GLUON_FILES ANISOTROPY
    EXTRA_SOURCES
    ../generic_pg/update.c
    ../generic_pg/update_u.c
    ../generic_pg/update_h.c
    ../generic/ranmom.c
)

add_milc_executable(su3_hmc control.c
    DEFINES HMC_ALGORITHM ONLY_GLUON_FILES
    EXTRA_SOURCES
    ../generic_pg/update.c
    ../generic_pg/update_u.c
    ../generic_pg/update_h.c
    ../generic/ranmom.c
    ../generic_pg/d_action.c
)

add_milc_executable(su3_hmc_a control.c
    DEFINES HMC_ALGORITHM ONLY_GLUON_FILES ANISOTROPY
    EXTRA_SOURCES
    ../generic_pg/update.c
    ../generic_pg/update_u.c
    ../generic_pg/update_h.c
    ../generic/ranmom.c
    ../generic_pg/d_action.c
)

add_milc_executable(su3_ora control.c
    DEFINES ORA_ALGORITHM ONLY_GLUON_FILES
    EXTRA_SOURCES
    ../generic_pg/update_ora.c
    ../generic_pg/dsdu_qhb.c
    ../generic_pg/relax.c
    ../generic_pg/monte.c
    ../generic/gaugefix2.c
)

add_milc_executable(su3_ora_a control.c
    DEFINES ORA_ALGORITHM ONLY_GLUON_FILES ANISOTROPY
    EXTRA_SOURCES
    ../generic_pg/update_ora.c
    ../generic_pg/dsdu_qhb.c
    ../generic_pg/relax.c
    ../generic_pg/monte.c
    ../generic/gaugefix2.c
)

add_milc_executable(su3_ora_glue control.c
    DEFINES ORA_ALGORITHM GBALL ONLY_GLUON_FILES
    EXTRA_SOURCES
    ../generic_pg/update_ora.c
    ../generic_pg/dsdu_qhb.c
    ../generic_pg/relax.c
    ../generic_pg/monte.c
    ../generic/gaugefix2.c
    ../generic/glueball_op.c
    ../generic/path_product.c
)
