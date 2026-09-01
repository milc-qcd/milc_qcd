add_interface_library(my_objects
    gauge_info.c
    staple.c
    integrate.c
    fmunu.c
    setup.c
)

add_interface_library(g_objects
    ../generic/check_unitarity.c
    ../generic/d_plaq4.c
    ../generic/field_strength.c
    ../generic/field_translation.c
    ../generic/field_utilities.c
    ../generic/gauge_utilities.c
    ../generic/gaugefix2.c
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
    ../generic/show_generic_opts.c
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

add_milc_executable(wilson_flow control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_LUSCHER
)

add_milc_executable(wilson_flow_a control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_LUSCHER ANISOTROPY
)

add_milc_executable(wilson_flow_cf3 control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_CF3
)

add_milc_executable(wilson_flow_cf3_a control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_CF3 ANISOTROPY
)

add_milc_executable(wilson_flow_ck control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_CK
)

add_milc_executable(wilson_flow_ck_a control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_CK ANISOTROPY
)

add_milc_executable(wilson_flow_bbb control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_BBB
)

add_milc_executable(wilson_flow_bbb_a control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_BBB ANISOTROPY
)

add_milc_executable(wilson_flow_rkmk3 control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_RKMK3
)

add_milc_executable(wilson_flow_rkmk3_a control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_RKMK3 ANISOTROPY
)

add_milc_executable(wilson_flow_rkmk4 control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_RKMK4
)

add_milc_executable(wilson_flow_rkmk4_a control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_RKMK4 ANISOTROPY
)

add_milc_executable(wilson_flow_rkmk5 control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_RKMK5
)

add_milc_executable(wilson_flow_rkmk5_a control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_RKMK5 ANISOTROPY
)

add_milc_executable(wilson_flow_rkmk8 control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_RKMK8
)

add_milc_executable(wilson_flow_rkmk8_a control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_RKMK8 ANISOTROPY
)

add_milc_executable(wilson_flow_adpt control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_ADAPT_LUSCHER
)

add_milc_executable(wilson_flow_adpt_a control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_ADAPT_LUSCHER ANISOTROPY
)

add_milc_executable(wilson_flow_adpt_cf3 control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_ADAPT_CF3
)

add_milc_executable(wilson_flow_adpt_cf3_a control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_ADAPT_CF3 ANISOTROPY
)

add_milc_executable(wilson_flow_adpt_bs control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_ADAPT_BS
)

add_milc_executable(wilson_flow_adpt_bs_a control.c
    DEFINES GF_INTEGRATOR=INTEGRATOR_ADAPT_BS ANISOTROPY
)
