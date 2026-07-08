add_library(misc_objs OBJECT
    crc32.c
)

set(complex_sources
    cadd.c cdiv.c ce_itheta.c cexp.c clog.c cmplx.c cmul.c
    conjg.c csqrt.c csub.c dcadd.c dcdiv.c dce_itheta.c
    dcexp.c dclog.c dcmplx.c dcmul.c dconjg.c dcsqrt.c dcsub.c
)

set(su3_sources
    addmat.c addvec.c cmp_ahmat.c cs_m_a_vec.c cs_m_a_mat.c
    cs_m_s_vec.c cs_m_vec.c det_su3.c clear_mat.c dumpmat.c
    dumptensor4.c dumpvec.c clearvec.c m_amatvec_s.c
    m_amatvec.c m_amatvec_ns.c
    m_mat_an.c m_mat_na.c m_mat_nn.c m_matvec.c m_matvec_ns.c
    m_matvec_s.c make_ahmat.c rand_ahmat.c realtr.c complextr.c
    s_m_a_hwvec_proj.c s_m_a_mat.c s_m_a_vec.c s_m_s_mat.c
    s_m_s_vec.c s_m_sum_vec.c s_m_vec.c s_m_mat.c cs_m_mat.c
    cs_m_s_mat.c su3_adjoint.c su3_dot.c su3_rdot.c su3_proj.c
    su3mat_copy.c su3vec_copy.c submat.c subvec.c trace_su3.c
    uncmp_ahmat.c msq_su3vec.c sub4vecs.c m_amv_4dir.c
    m_amv_4vec.c m_mv_s_4dir.c flush_to_zero.c
    l_su2_hit_n.c r_su2_hit_a.c dumpsu2.c
    m_su2_mat_vec_n.c m_su2_mat_vec_a.c
    gaussrand.c z2rand.c byterevn.c prefetch.c eigen_su3_UdU.c
)

set(wilson_sources
    wp_shrink.c wp_grow.c wp_grow_a.c dump_wvec.c
    clear_wvec.c su3_proj_w.c copy_wvec.c add_wvec.c sub_wvec.c
    s_m_wvec.c s_m_hwvec.c msq_wvec.c wvec_dot.c wvec2_dot.c
    wvec_rdot.c s_m_a_wvec.c s_m_atm_wvec.c mb_gamma.c mb_gamma_l.c
    mb_gamma_r.c mswvb_gamma_l.c mswvb_gamma_r.c
    cs_m_wvec.c cs_m_a_wvec.c cs_m_a_wvec2.c
    m_mat_wvec.c m_mat_hwvec.c m_amat_wvec.c m_amat_hwvec.c
    grow4wvecs.c wp_shrink4.c
)

add_library(qcd INTERFACE)
if(PRECISION EQUAL 1)
    add_library(complex.1 STATIC ${complex_sources} $<TARGET_OBJECTS:misc_objs>)
    target_compile_definitions(complex.1 PRIVATE MILC_PRECISION=1)
    add_library(su3.1 STATIC ${su3_sources} ${wilson_sources})
    target_compile_definitions(su3.1 PRIVATE MILC_PRECISION=1)
    target_link_libraries(qcd INTERFACE su3.1 complex.1)
elseif(PRECISION EQUAL 2)
    add_library(complex.2 STATIC ${complex_sources} $<TARGET_OBJECTS:misc_objs>)
    target_compile_definitions(complex.2 PRIVATE MILC_PRECISION=2)
    add_library(su3.2 STATIC ${su3_sources} ${wilson_sources})
    target_compile_definitions(su3.2 PRIVATE MILC_PRECISION=2)
    target_link_libraries(qcd INTERFACE su3.2 complex.2)
endif()
