# Source groups and targets follow Make_template.
# Legacy machine-dependent I/O uses the current serial/MPI communication
# and ANSI I/O interfaces; old source names are retained below.

add_interface_library(my_objects
)

add_interface_library(g_objects
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

add_milc_executable(dense_su3 control_dense.c
    DEFINES ORA_ALGORITHM
    EXTRA_SOURCES
    setup.c
    ../generic/ploop3.c
    ../generic/ranstuff.c
    # TODO: Restore the source referenced by Make_template.
    layout_hyper.c
    ../generic/io_lat4.c
    ../generic/plaquette4.c
    ../generic/check_unitarity.c
    ../generic/make_lattice.c
    ../generic/reunitarize2.c
    dsdu_qhb.c
    relax_space.c
    metropolis_dense.c
    update_dense.c
    density.c
    heatbath_space.c
    ../generic/io_helpers.c
    gauge_info.c
    ../generic_ks/rephase.c
)

add_milc_executable(dense_eff_su3 control_dense.c
    DEFINES ORA_ALGORITHM AXIAL_GAUGE EFFECTIVE
    EXTRA_SOURCES
    setup.c
    ../generic/ploop3.c
    ../generic/ranstuff.c
    # TODO: Restore the source referenced by Make_template.
    layout_hyper.c
    ../generic/io_lat4.c
    ../generic/plaquette4.c
    ../generic/check_unitarity.c
    ../generic/make_lattice.c
    ../generic/reunitarize2.c
    dsdu_qhb.c
    relax_space.c
    metropolis_dense.c
    update_dense.c
    density.c
    heatbath_space.c
    ../generic/io_helpers.c
    gauge_info.c
    ../generic_ks/rephase.c
    # TODO: Restore the source referenced by Make_template.
    density_effective.c
    # TODO: Restore the source referenced by Make_template.
    external_field.c
)

add_milc_executable(halfdense_su3 control_dense.c
    DEFINES ORA_ALGORITHM
    EXTRA_SOURCES
    setup.c
    ../generic/ploop3.c
    ../generic/ranstuff.c
    # TODO: Restore the source referenced by Make_template.
    layout_hyper.c
    ../generic/plaquette4.c
    ../generic/check_unitarity.c
    ../generic/make_lattice.c
    ../generic/reunitarize2.c
    dsdu_qhb.c
    relax_space.c
    metropolis_halfdense.c
    update_dense.c
    density_half.c
    heatbath_space.c
    ../generic/io_helpers.c
    gauge_info.c
    ../generic_ks/rephase.c
)

add_milc_executable(dense_axial_su3 control_dense.c
    DEFINES ORA_ALGORITHM AXIAL_GAUGE
    EXTRA_SOURCES
    setup.c
    ../generic/ploop3.c
    ../generic/ranstuff.c
    # TODO: Restore the source referenced by Make_template.
    layout_hyper.c
    ../generic/io_lat4.c
    ../generic/plaquette4.c
    ../generic/check_unitarity.c
    ../generic/make_lattice.c
    ../generic/reunitarize2.c
    dsdu_qhb.c
    relax_space.c
    metropolis_dense.c
    update_dense.c
    density.c
    heatbath_space.c
    ../generic/io_helpers.c
    gauge_info.c
    ../generic_ks/rephase.c
)

add_milc_executable(halfdense_axial_su3 control_dense.c
    DEFINES ORA_ALGORITHM AXIAL_GAUGE
    EXTRA_SOURCES
    setup.c
    ../generic/ploop3.c
    ../generic/ranstuff.c
    # TODO: Restore the source referenced by Make_template.
    layout_hyper.c
    ../generic/plaquette4.c
    ../generic/check_unitarity.c
    ../generic/make_lattice.c
    ../generic/reunitarize2.c
    dsdu_qhb.c
    relax_space.c
    metropolis_halfdense.c
    update_dense.c
    density_half.c
    heatbath_space.c
    ../generic/io_helpers.c
    gauge_info.c
    ../generic_ks/rephase.c
)

add_milc_executable(dense_su3_lpbp control_dense.c
    DEFINES LIGHT_PBP ORA_ALGORITHM
    EXTRA_SOURCES
    setup.c
    ../generic/ploop3.c
    ../generic/ranstuff.c
    # TODO: Restore the source referenced by Make_template.
    layout_hyper.c
    ../generic/io_lat4.c
    ../generic/plaquette4.c
    ../generic/check_unitarity.c
    ../generic/make_lattice.c
    ../generic/reunitarize2.c
    dsdu_qhb.c
    relax_space.c
    metropolis_dense.c
    update_dense.c
    density.c
    heatbath_space.c
    ../generic/io_helpers.c
    gauge_info.c
    ../generic_ks/rephase.c
    ../generic_ks/d_congrad5.c
    ../generic_ks/grsource.c
    f_measure.c
)

