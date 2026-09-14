# Source groups and targets follow Make_template.

add_interface_library(my_objects
    setup.c
    gauge_info.c
)

add_interface_library(g_objects
    ../generic/layout_hyper_prime.c
)

add_interface_library(objects
    OBJECTS
    my_objects
    g_objects
    fft_objects
    gauge_objects
    standard_objects
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

add_milc_executable(su3_gluon_prop_nogfix control.c
    DEFINES GLUON_PROP
    EXTRA_SOURCES
    gluon_prop.c
)

add_milc_executable(su3_gluon_prop_imp_nogfix control.c
    DEFINES GLUON_PROP
    EXTRA_SOURCES
    gluon_prop_imp.c
)

add_milc_executable(su3_gluon_prop control.c
    DEFINES GLUON_PROP GFIX
    EXTRA_SOURCES
    gluon_prop.c
    gaugefixfft.c
)

add_milc_executable(su3_gluon_prop_imp control.c
    DEFINES GLUON_PROP GFIX IMP_GFIX
    EXTRA_SOURCES
    gluon_prop_imp.c
    gaugefixfft.c
)

add_milc_executable(su3_ks_quark_prop control.c
    QUARK onelink_action
    DEFINES FN QUARK_PROP GFIX HAVE_KS
    EXTRA_OBJECTS fn_objects fermion_objects ks_io_objects
    EXTRA_SOURCES
    ksprop_info.c
    ../generic/ploop3.c
    quark_prop.c
    gaugefixfft.c
)

add_milc_executable(su3_fn_quark_prop control.c
    QUARK naik_action
    DEFINES FN QUARK_PROP GFIX NAIK HAVE_KS
    EXTRA_OBJECTS fn_objects fermion_objects ks_io_objects
    EXTRA_SOURCES
    ksprop_info.c
    ../generic/ploop3.c
    quark_prop.c
    gaugefixfft.c
)

add_milc_executable(su3_p4_quark_prop control.c
    QUARK p4_action
    DEFINES EO QUARK_PROP GFIX HAVE_KS
    EXTRA_OBJECTS eo_objects multi_inv_eo fermion_objects ks_io_objects
    EXTRA_SOURCES
    ksprop_info.c
    ../generic/ploop3.c
    quark_prop.c
    gaugefixfft.c
)

add_milc_executable(su3_asqtad_quark_prop control.c
    QUARK asqtad_action
    DEFINES FN QUARK_PROP GFIX NAIK HAVE_KS
    EXTRA_OBJECTS fn_objects fermion_objects ks_io_objects
    EXTRA_SOURCES
    ksprop_info.c
    ../generic/ploop3.c
    quark_prop.c
    gaugefixfft.c
)

add_milc_executable(su3_asqtad_renorm control.c
    QUARK asqtad_action
    DEFINES FN QUARK_PROP QUARK_RENORM GFIX NAIK HAVE_KS
    EXTRA_OBJECTS fn_objects fermion_objects ks_io_objects
    EXTRA_SOURCES
    ksprop_info.c
    ../generic/ploop3.c
    quark_renorm.c
    gaugefixfft.c
)

