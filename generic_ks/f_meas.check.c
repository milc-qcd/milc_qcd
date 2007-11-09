/**************** f_meas.c ***************************************/
/* MIMD version 7 */
/* UH 11/01 Include measurements for quark number susceptibilities */
/*          Note: this does not work for p4-action! */
/* UH 11/1/01 write complex stochastic estimators */
/* TB 10/01 Include measurements of dM/du0 for EOS */
/* CD 7/14/01 allow for multiple stochastic estimators NPBP_REPS */
/* DT 12/97 */
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.
*/

/* Measure fermionic observables:
    psi-bar-psi (separately on even and odd sites)
    fermion action

    This routine uses g_rand, phi, xxx, and other vectors used by
    the matrix inversion.
*/

#include "generic_ks_includes.h"	/* definitions files and prototypes */

void f_meas_imp( field_offset phi_off, field_offset xxx_off, Real mass,
		 ferm_links_t *fn, ferm_links_t *fn_dmdu0 ){
    Real r_psi_bar_psi_even, i_psi_bar_psi_even;
    Real  r_psi_bar_psi_odd, i_psi_bar_psi_odd;
    Real r_ferm_action;
    /* local variables for accumulators */
    register int i;
    register site *st;
    double rfaction;
    double_complex pbp_e, pbp_o;
#ifdef NPBP_REPS
    double pbp_pbp;
#endif
    complex cc;

#ifdef DM_DU0
    double r_pb_dMdu_p_even, r_pb_dMdu_p_odd;
    /* check */ double r_gb_dMdu_g_e, r_gb_dMdu_g_o;
    /* check */ double r_gb_M_g_e, r_gb_M_g_o;
#endif

#ifdef CHEM_POT
#ifndef FN		/* FN is assumed for quark number susc. */
BOMB THE COMPILE
#endif
#ifndef NPBP_REPS	/* Need multiple repetitions for susceptibilities! */
BOMB THE COMPILE
#endif
    msg_tag *tag0, *tag1, *tag2, *tag3;
    double_complex pb_dMdmu_p_e, pb_dMdmu_p_o;
    double_complex pb_d2Mdmu2_p_e, pb_d2Mdmu2_p_o;
    Real r_pb_dMdmu_p_e, i_pb_dMdmu_p_e;
    Real r_pb_dMdmu_p_o, i_pb_dMdmu_p_o;
    Real r_pb_d2Mdmu2_p_e, i_pb_d2Mdmu2_p_e;
    Real r_pb_d2Mdmu2_p_o, i_pb_d2Mdmu2_p_o;
    double MidM_MidM;
#endif

#ifdef NPBP_REPS
    int npbp_reps = npbp_reps_in;  /* Number of repetitions of stochastic
                                   estimate */
    int prec = prec_pbp;  /* Precision of the inversion */
#else
    int npbp_reps = 1;   /* Default values */
    int prec = PRECISION;
#endif
    int jpbp_reps;
#endif

    for(jpbp_reps = 0; jpbp_reps < npbp_reps; jpbp_reps++){
      rfaction = (double)0.0;
      pbp_e = pbp_o = dcmplx((double)0.0,(double)0.0);
      
      /* Make random source, and do inversion */
      /* generate g_rand random; phi_off = M g_rand */
      grsource_imp( phi_off, mass, EVENANDODD );
      /* phi_off = M g_rand (still) */
      /* xxx_off = M^{-1} g_rand */
      //      clear_latvec(xxx_off,EVENANDODD);
      mat_invert_uml( F_OFFSET(g_rand), xxx_off, phi_off, mass, prec, fn );
      
#ifdef DM_DU0
      r_pb_dMdu_p_even = r_pb_dMdu_p_odd = (double)0.0;
      /* dMdu_x = dM/du0 M^{-1} g_rand */
      ddslash_fn_du0_site( xxx_off, F_OFFSET(dMdu_x), EVENANDODD, 
			   fn_dmdu0 );
      /* check */ r_gb_dMdu_g_e = r_gb_dMdu_g_o = 0;
      /* check */ r_gb_M_g_e = r_gb_M_g_o = 0;
      /* check */ ddslash_fn_du0_site (F_OFFSET(g_rand), F_OFFSET(dM_check), EVENANDODD, fn_dmdu0 );
#endif

#ifdef CHEM_POT
      pb_dMdmu_p_e = pb_dMdmu_p_o = dcmplx((double)0.0,(double)0.0);
      pb_d2Mdmu2_p_e = pb_d2Mdmu2_p_o = dcmplx((double)0.0,(double)0.0);

      /* Start gathers from positive t-direction */
      tag0 = start_gather_site( xxx_off, sizeof(su3_vector), TUP,
	EVENANDODD, gen_pt[0] );
      tag1 = start_gather_site( xxx_off, sizeof(su3_vector), T3UP,
	EVENANDODD, gen_pt[1] );

      FORALLSITES(i,st){
	mult_adj_su3_mat_vec( &(t_fatlink[4*i+TUP]),
	    (su3_vector *)F_PT(st,xxx_off), &(st->tempvec[TUP]) );
	mult_adj_su3_mat_vec( &(t_longlink[4*i+TUP]),
	    (su3_vector *)F_PT(st,xxx_off), &(st->templongvec[TUP]) );
      }

      /* Start gathers from negative t-direction */
      tag2 = start_gather_site( F_OFFSET(tempvec[TUP]), sizeof(su3_vector),
	OPP_DIR(TUP), EVENANDODD, gen_pt[2] );
      tag3 = start_gather_site( F_OFFSET(templongvec[TUP]), sizeof(su3_vector),
	OPP_3_DIR(T3UP), EVENANDODD, gen_pt[3] );

      /* Wait gathers from positive t-direction and multiply by matrix */
      wait_gather(tag0);
      wait_gather(tag1);

      FORALLSITES(i,st){
	mult_su3_mat_vec( &(t_fatlink[4*i+TUP]),
	    (su3_vector *)gen_pt[0][i], &(st->tempvec[0]) );
	mult_su3_mat_vec( &(t_longlink[4*i+TUP]),
	    (su3_vector *)gen_pt[1][i], &(st->templongvec[0]) );
      }

      /* Wait gathers from negative t-direction */
      wait_gather(tag2);
      wait_gather(tag3);
#endif

      /* fermion action = phi.xxx */
      /* psi-bar-psi on even sites = g_rand.xxx */
      FOREVENSITES(i,st){
	cc = su3_dot( (su3_vector *)F_PT(st,phi_off),
		      (su3_vector *)F_PT(st,xxx_off) );
	rfaction += cc.real;
	cc = su3_dot( &(st->g_rand), (su3_vector *)F_PT(st,xxx_off) );
	CSUM(pbp_e, cc);
#ifdef DM_DU0
	/* r_pb_dMdu_p_even = g_rand * dM/du0 M^{-1} g_rand |even*/
	cc = su3_dot( &(st->g_rand), &(st->dMdu_x) );
	r_pb_dMdu_p_even += cc.real;
	/* check */ cc = su3_dot( &(st->g_rand), &(st->dM_check) );
	/* check */ r_gb_dMdu_g_e += cc.real;
	/* check */ cc = su3_dot( &(st->g_rand), (su3_vector *)F_PT(st,phi_off) );
	/* check */ r_gb_M_g_e += cc.real;
#endif

#ifdef CHEM_POT
	cc = su3_dot( &(st->g_rand), &(st->tempvec[0]) );
	CSUM(pb_dMdmu_p_e, cc);
	CSUM(pb_d2Mdmu2_p_e, cc);
	cc = su3_dot( &(st->g_rand), (su3_vector *)gen_pt[2][i] );
	CSUM(pb_dMdmu_p_e, cc);
	CSUB(pb_d2Mdmu2_p_e, cc, pb_d2Mdmu2_p_e);
	cc = su3_dot( &(st->g_rand), &(st->templongvec[0]) );
	CMULREAL(cc, 3.0, cc);
	CSUM(pb_dMdmu_p_e, cc);
	CMULREAL(cc, 3.0, cc);
	CSUM(pb_d2Mdmu2_p_e, cc);
	cc = su3_dot( &(st->g_rand), (su3_vector *)gen_pt[3][i] );
	CMULREAL(cc, 3.0, cc);
	CSUM(pb_dMdmu_p_e, cc);
	CMULREAL(cc, 3.0, cc);
	CSUB(pb_d2Mdmu2_p_e, cc, pb_d2Mdmu2_p_e);
	add_su3_vector( &(st->tempvec[0]), (su3_vector *)gen_pt[2][i],
			&(st->tempvec[0]) );
	add_su3_vector( &(st->templongvec[0]), (su3_vector *)gen_pt[3][i],
			&(st->templongvec[0]) );
	scalar_mult_add_su3_vector( &(st->tempvec[0]), &(st->templongvec[0]),
				    3.0, &(st->dM_M_inv) );
#endif
      }

      /* psi-bar-psi on odd sites */
      FORODDSITES(i,st){
	cc = su3_dot( &(st->g_rand), (su3_vector *)F_PT(st,xxx_off) );
	CSUM(pbp_o, cc);
#ifdef DM_DU0
	/* r_pb_dMdu_p_odd = g_rand * dM/du0 M^{-1} g_rand |odd*/
	cc = su3_dot( &(st->g_rand), &(st->dMdu_x) );
	r_pb_dMdu_p_odd += cc.real;
	/* check */ cc = su3_dot( &(st->g_rand), &(st->dM_check) );
	/* check */ r_gb_dMdu_g_o += cc.real;
	/* check */ cc = su3_dot( &(st->g_rand), (su3_vector *)F_PT(st,phi_off) );
	/* check */ r_gb_M_g_o += cc.real;
#endif

#ifdef CHEM_POT
	cc = su3_dot( &(st->g_rand), &(st->tempvec[0]) );
	CSUM(pb_dMdmu_p_o, cc);
	CSUM(pb_d2Mdmu2_p_o, cc);
	cc = su3_dot( &(st->g_rand), (su3_vector *)gen_pt[2][i] );
	CSUM(pb_dMdmu_p_o, cc);
	CSUB(pb_d2Mdmu2_p_o, cc, pb_d2Mdmu2_p_o);
	cc = su3_dot( &(st->g_rand), &(st->templongvec[0]) );
	CMULREAL(cc, 3.0, cc);
	CSUM(pb_dMdmu_p_o, cc);
	CMULREAL(cc, 3.0, cc);
	CSUM(pb_d2Mdmu2_p_o, cc);
	cc = su3_dot( &(st->g_rand), (su3_vector *)gen_pt[3][i] );
	CMULREAL(cc, 3.0, cc);
	CSUM(pb_dMdmu_p_o, cc);
	CMULREAL(cc, 3.0, cc);
	CSUB(pb_d2Mdmu2_p_o, cc, pb_d2Mdmu2_p_o);
	add_su3_vector( &(st->tempvec[0]), (su3_vector *)gen_pt[2][i],
			&(st->tempvec[0]) );
	add_su3_vector( &(st->templongvec[0]), (su3_vector *)gen_pt[3][i],
			&(st->templongvec[0]) );
	scalar_mult_add_su3_vector( &(st->tempvec[0]), &(st->templongvec[0]),
				    3.0, &(st->dM_M_inv) );
#endif
      }
      
      g_dcomplexsum( &pbp_o );
      g_dcomplexsum( &pbp_e );
      g_doublesum( &rfaction );
      
#ifdef DM_DU0
      g_doublesum( &r_pb_dMdu_p_even );
      g_doublesum( &r_pb_dMdu_p_odd );
      r_pb_dMdu_p_even *= (2.0/(double)volume);
      r_pb_dMdu_p_odd *= (2.0/(double)volume);
      node0_printf("PB_DMDU_P: mass %e  %e  %e ( %d of %d )\n", mass,
		   r_pb_dMdu_p_even, r_pb_dMdu_p_odd, jpbp_reps+1, npbp_reps);
      /* check */ g_doublesum( &r_gb_dMdu_g_e );
      /* check */ g_doublesum( &r_gb_dMdu_g_o );
      /* check */ r_gb_dMdu_g_e *= (2.0/(double)volume);
      /* check */ r_gb_dMdu_g_o *= (2.0/(double)volume);
      /* check */ node0_printf("G_DMDU_G: mass %e  %e  %e ( %d of %d )\n", mass, r_gb_dMdu_g_e, r_gb_dMdu_g_o, jpbp_reps+1, npbp_reps);
      /* check */ g_doublesum( &r_gb_M_g_e );
      /* check */ g_doublesum( &r_gb_M_g_o );
      /* check */ r_gb_M_g_e *= (2.0/(double)volume);
      /* check */ r_gb_M_g_o *= (2.0/(double)volume);
      /* check */ node0_printf("G_M_G: mass %e  %e  %e ( %d of %d )\n", mass, r_gb_M_g_e, r_gb_M_g_o, jpbp_reps+1, npbp_reps);

#endif

      r_psi_bar_psi_odd =  pbp_o.real*(2.0/(double)volume) ;
      i_psi_bar_psi_odd =  pbp_o.imag*(2.0/(double)volume) ;
      r_psi_bar_psi_even =  pbp_e.real*(2.0/(double)volume) ;
      i_psi_bar_psi_even =  pbp_e.imag*(2.0/(double)volume) ;
      r_ferm_action =  rfaction*(1.0/(double)volume) ;
      node0_printf("PBP: mass %e     %e  %e  %e  %e ( %d of %d )\n", mass,
		   r_psi_bar_psi_even, r_psi_bar_psi_odd,
		   i_psi_bar_psi_even, i_psi_bar_psi_odd,
		   jpbp_reps+1, npbp_reps);
      node0_printf("FACTION: mass = %e,  %e ( %d of %d )\n", mass,
		   r_ferm_action, jpbp_reps+1, npbp_reps);

#ifdef CHEM_POT
      /* free up the buffers */
      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
      cleanup_gather(tag3);

      g_dcomplexsum( &pb_dMdmu_p_e );
      g_dcomplexsum( &pb_dMdmu_p_o );
      g_dcomplexsum( &pb_d2Mdmu2_p_e );
      g_dcomplexsum( &pb_d2Mdmu2_p_o );

      r_pb_dMdmu_p_e =  pb_dMdmu_p_e.real*(2.0/(double)volume) ;
      i_pb_dMdmu_p_e =  pb_dMdmu_p_e.imag*(2.0/(double)volume) ;
      r_pb_dMdmu_p_o =  pb_dMdmu_p_o.real*(2.0/(double)volume) ;
      i_pb_dMdmu_p_o =  pb_dMdmu_p_o.imag*(2.0/(double)volume) ;
      r_pb_d2Mdmu2_p_e =  pb_d2Mdmu2_p_e.real*(2.0/(double)volume) ;
      i_pb_d2Mdmu2_p_e =  pb_d2Mdmu2_p_e.imag*(2.0/(double)volume) ;
      r_pb_d2Mdmu2_p_o =  pb_d2Mdmu2_p_o.real*(2.0/(double)volume) ;
      i_pb_d2Mdmu2_p_o =  pb_d2Mdmu2_p_o.imag*(2.0/(double)volume) ;
      node0_printf("PB_DMDMU_P: mass %e  %e  %e  %e  %e ( %d of %d )\n", mass,
		   r_pb_dMdmu_p_e, r_pb_dMdmu_p_o,
		   i_pb_dMdmu_p_e, i_pb_dMdmu_p_o,
		   jpbp_reps+1, npbp_reps);
      node0_printf("PB_D2MDMU2_P: mass %e  %e  %e  %e  %e ( %d of %d )\n", mass,
		   r_pb_d2Mdmu2_p_e, r_pb_d2Mdmu2_p_o,
		   i_pb_d2Mdmu2_p_e, i_pb_d2Mdmu2_p_o,
		   jpbp_reps+1, npbp_reps);
#endif

#ifdef NPBP_REPS
      pbp_pbp = (double)0.0;
      FORALLSITES(i,st){
	su3vec_copy( (su3_vector *)F_PT(st,xxx_off), &(st->M_inv) );
      }
      mat_invert_uml( F_OFFSET(M_inv), xxx_off, phi_off, mass, fn );
      FORALLSITES(i,st){
	cc = su3_dot( &(st->g_rand), (su3_vector *)F_PT(st,xxx_off) );
	pbp_pbp += cc.real;
      }
      g_doublesum( &pbp_pbp );
      pbp_pbp =  pbp_pbp*(1.0/(double)volume) ;
      node0_printf("TR_MM_INV: mass %e,  %e ( %d of %d )\n", mass,
		   pbp_pbp, jpbp_reps+1, npbp_reps);
#endif

#ifdef CHEM_POT
      mat_invert_uml( F_OFFSET(dM_M_inv), xxx_off, phi_off, mass, fn );

      /* Start gathers from positive t-direction */
      tag0 = start_gather_site( xxx_off, sizeof(su3_vector), TUP,
	EVENANDODD, gen_pt[0] );
      tag1 = start_gather_site( xxx_off, sizeof(su3_vector), T3UP,
	EVENANDODD, gen_pt[1] );

      FORALLSITES(i,st){
	mult_adj_su3_mat_vec( &(t_fatlink[4*i+TUP]),
	    (su3_vector *)F_PT(st,xxx_off), &(st->tempvec[TUP]) );
	mult_adj_su3_mat_vec( &(t_longlink[4*i+TUP]),
	    (su3_vector *)F_PT(st,xxx_off), &(st->templongvec[TUP]) );
      }

      /* Start gathers from negative t-direction */
      tag2 = start_gather_site( F_OFFSET(tempvec[TUP]), sizeof(su3_vector),
	OPP_DIR(TUP), EVENANDODD, gen_pt[2] );
      tag3 = start_gather_site( F_OFFSET(templongvec[TUP]), sizeof(su3_vector),
	OPP_3_DIR(T3UP), EVENANDODD, gen_pt[3] );

      /* Wait gathers from positive t-direction and multiply by matrix */
      wait_gather(tag0);
      wait_gather(tag1);

      FORALLSITES(i,st){
	mult_su3_mat_vec( &(t_fatlink[4*i+TUP]),
	    (su3_vector *)gen_pt[0][i], &(st->tempvec[0]) );
	mult_su3_mat_vec( &(t_longlink[4*i+TUP]),
	    (su3_vector *)gen_pt[1][i], &(st->templongvec[0]) );
      }

      /* Wait gathers from negative t-direction */
      wait_gather(tag2);
      wait_gather(tag3);

      MidM_MidM = (double)0.0;

      FORALLSITES(i,st){
	cc = su3_dot( &(st->g_rand), &(st->tempvec[0]) );
	MidM_MidM += cc.real;
	cc = su3_dot( &(st->g_rand), (su3_vector *)gen_pt[2][i] );
	MidM_MidM += cc.real;
	cc = su3_dot( &(st->g_rand), &(st->templongvec[0]) );
	MidM_MidM += 3.0 * cc.real;
	cc = su3_dot( &(st->g_rand), (su3_vector *)gen_pt[3][i] );
	MidM_MidM += 3.0 * cc.real;
      }

      g_doublesum( &MidM_MidM );
      MidM_MidM =  MidM_MidM*(1.0/(double)volume) ;
      node0_printf("TR_MidM_MidM: mass %e,  %e ( %d of %d )\n", mass,
		   MidM_MidM, jpbp_reps+1, npbp_reps);

      /* free up the buffers */
      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
      cleanup_gather(tag3);
#endif
    }
}

