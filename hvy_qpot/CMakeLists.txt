add_interface_library(layout
    ../generic/layout_hyper_prime.c
)

add_interface_library(my_objects
    gauge_info.c
    setup.c
)

add_interface_library(g_ks_objects
    ../generic_ks/fermion_links_hyp.c
)

add_interface_library(add_objects
    ../generic/hvy_pot.c
    OBJECTS
    machine_dep_io
    commpkg
    g_ks_objects
)

add_interface_library(objects
    OBJECTS
    my_objects
    gauge_objects
    layout
    standard_objects
    add_objects
)

add_interface_library(libraries
    OBJECTS
    qcd
    scidac
    add
    m
)

add_milc_executable(su3_hqp_ape control.c
    EXTRA_SOURCES
    ../generic/ax_gauge.c
    gball_simp.c
    hybrid_loop1.c
    smearing.c
    w_loop1.c
    w_loop2.c

)

add_milc_executable(su3_hqp_hyp control.c
    DEFINES HYP_3D_SMEARING
    EXTRA_SOURCES
    ../generic/ax_gauge.c
    gball_simp.c
    hybrid_loop1.c
    smearing_hyp.c
    w_loop1.c
    w_loop2.c
)

add_milc_executable(su3_hybrid_hqp_ape control.c
    DEFINES HYBRIDS_MEASURE
    EXTRA_SOURCES
    ../generic/ax_gauge.c
    gball_simp.c
    hybrid_loop1.c
    smearing.c
    w_loop1.c
    w_loop2.c
)

add_milc_executable(su3_hybrid_hqp_hyp control.c
    DEFINES HYBRIDS_MEASURE HYP_3D_SMEARING
    EXTRA_SOURCES
    ../generic/ax_gauge.c
    gball_simp.c
    hybrid_loop1.c
    smearing_hyp.c
    w_loop1.c
    w_loop2.c
)

add_milc_executable(su3_hqp_coulomb control.c
    DEFINES COULOMB
)

add_milc_executable(su3_hqp_coulomb_plcor control.c
    DEFINES COULOMB PLOOPCOR_MEAS
)
