/**************** f_meas.c ***************************************/
/* MIMD version 7 */
/* CD 12/10 Added multimass version and consolidated code */ 
/* UH 11/01 Include measurements for quark number susceptibilities */
/*          Note: this does not work for p4-action! */
/* UH 11/1/01 write complex stochastic estimators */
/* TB 10/01 Include measurements of dM/du0 for EOS */
/* CD 7/14/01 allow for multiple stochastic estimators NPBP_REPS -
   This became the default 09/27/09 CD */
/* DT 12/97 */
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.
*/

/* Measure fermionic observables:
    psi-bar-psi (separately on even and odd sites)
    fermion action

    With DM_DU0:
    pb_dMdu_p (derivative of pbp with respect to u0)

    With CHEM_POT:
    pb_dMdmu_p 
    pb_d2Mdmu2_p 
    Tr_MidM_MidM

    Entry points

    f_meas_imp_field
    f_meas_imp_multi
    f_meas_imp (deprecated)
*/

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/fn_links.h"

#if FERM_ACTION == HISQ & defined(DM_DU0)
BOMB THE COMPILATION No dM_du0 support for HISQ actions
#endif

static Real *create_real_array(int n){
  Real *a;
  int i;

  a = (Real *)malloc(n*sizeof(Real));
  if(a == NULL){
    printf("f_meas: No room for array\n");
    terminate(1);
  }
  for(i = 0; i < n; i++)a[i] = 0.;
  return a;
}

static void destroy_real_array(Real *a){
  if(a == NULL)return;
  free(a);
}

static su3_vector **create_su3_vector_array(int n){
  su3_vector **a;
  int i;

  a = (su3_vector **)malloc(n*sizeof(su3_vector *));
  if(a == NULL){
    printf("f_meas: No room for array\n");
    terminate(1);
  }
  for(i = 0; i < n; i++) a[i] = create_v_field();
  return a;
}

static void destroy_su3_vector_array(su3_vector **a, int n){
  int i;

  if(a == NULL)return;
  for(i = 0; i < n; i++)
    if(a[i] != NULL)
      destroy_v_field(a[i]);
}

#ifdef CHEM_POT

/* Act on the color vector field A with the time-like contribution to
   Dslash, but weight the long-link term by an additional factor of
   wtlong and multiply the entire backward contribution by sign.
   Result in dM_A. */

static void chem_pot_tshift(imp_ferm_links_t *fn, su3_vector *dM_A,
			    su3_vector *A, 
			    double wtlong, double sign){

  clear_v_field(dM_A);

  /* Parallel transport A from TUP via Dslash with
     weight wtlong for the long link contribution */
  dslash_fn_dir(A, dM_A, EVENANDODD, fn, TUP, +1, 1., wtlong );

  /* Add the contribution from the backward links with sign */
  dslash_fn_dir(A, dM_A, EVENANDODD, fn, TUP, -1, sign, sign*wtlong);
}  


static void chem_pot_print1(double_complex pb_dMdmu_p_e, double_complex pb_dMdmu_p_o,
			    double_complex pb_d2Mdmu2_p_e, double_complex pb_d2Mdmu2_p_o,
			    Real mass, int jpbp_reps, int npbp_reps)
{

  Real r_pb_dMdmu_p_e = 0.0;
  Real i_pb_dMdmu_p_e = 0.0;
  Real r_pb_dMdmu_p_o  = 0.0;
  Real i_pb_dMdmu_p_o = 0.0;
  Real r_pb_d2Mdmu2_p_e = 0.0;
  Real i_pb_d2Mdmu2_p_e = 0.0;
  Real r_pb_d2Mdmu2_p_o = 0.0;
  Real i_pb_d2Mdmu2_p_o = 0.0;

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

}

  
static void chem_pot_print2(double MidM_MidM, int jpbp_reps, int npbp_reps, 
			    Real mass){

  g_doublesum( &MidM_MidM );
  MidM_MidM =  MidM_MidM*(1.0/(double)volume) ;
  node0_printf("TR_MidM_MidM: mass %e,  %e ( %d of %d )\n", mass,
	       MidM_MidM, jpbp_reps+1, npbp_reps);
}

#endif

void f_meas_imp_field( int npbp_reps, quark_invert_control *qic, Real mass,
		       int naik_term_epsilon_index, fermion_links_t *fl){

  imp_ferm_links_t* fn = get_fm_links(fl)[naik_term_epsilon_index];

#ifdef DM_DU0
  imp_ferm_links_t* fn_du0 = get_fm_du0_links(fl)[naik_term_epsilon_index];
#endif

#if FERM_ACTION == HISQ & defined(DM_DEPS)
  imp_ferm_links_t *fn_deps = get_fn_deps_links(fl);
#endif

    Real r_psi_bar_psi_even, i_psi_bar_psi_even;
    Real  r_psi_bar_psi_odd, i_psi_bar_psi_odd;
    Real r_ferm_action;
    /* local variables for accumulators */
    register int i;
    double rfaction;
    double_complex pbp_e, pbp_o;
    complex cc;

    int jpbp_reps;
    su3_vector *gr = NULL;
    su3_vector *M_gr = NULL;
    su3_vector *M_inv_gr = NULL;

#ifdef DM_DU0
    double r_pb_dMdu_p_even, r_pb_dMdu_p_odd;
    su3_vector *dMdu_x = NULL;
#endif

#if FERM_ACTION == HISQ & defined(DM_DEPS)
    double r_pb_dMdeps_p_even, r_pb_dMdeps_p_odd;
    su3_vector *dMdeps_x = NULL;
#endif

#ifdef CHEM_POT
    double_complex pb_dMdmu_p_e, pb_dMdmu_p_o;
    double_complex pb_d2Mdmu2_p_e, pb_d2Mdmu2_p_o;
    double MidM_MidM;
    su3_vector *dM_M_inv_gr = NULL;
    su3_vector *d2M_M_inv_gr = NULL;
    su3_vector *M_inv_dM_M_inv_gr = NULL;
    su3_vector *dM_M_inv_dM_M_inv_gr = NULL;
#endif

    /* Loop over random sources */
    for(jpbp_reps = 0; jpbp_reps < npbp_reps; jpbp_reps++){

      rfaction = (double)0.0;
      pbp_e = pbp_o = dcmplx((double)0.0,(double)0.0);
      
      /* Make random source, and do inversion */
      /* generate gr random; M_gr = M gr */
      gr = create_v_field();
#ifndef Z2RSOURCE
      grsource_plain_field( gr, EVENANDODD );
#else
      z2rsource_plain_field( gr, mass, EVENANDODD );
#endif
      /* The following operation is done in the prevailing
	 precision.  The algorithm needs to be fixed! */
      M_gr = create_v_field();
      ks_dirac_adj_op( gr, M_gr, mass, EVENANDODD, fn );

      /* M_inv_gr = M^{-1} gr */

      M_inv_gr = create_v_field();
      mat_invert_uml_field( gr, M_inv_gr, qic, mass, fn );
      
#ifdef DM_DU0
      r_pb_dMdu_p_even = r_pb_dMdu_p_odd = (double)0.0;
      /* dMdu_x = dM/du0 M^{-1} gr */
      dMdu_x = create_v_field();
      dslash_fn_field( M_inv_gr, dMdu_x, EVENANDODD, fn_du0 );
#endif

#if FERM_ACTION == HISQ & defined(DM_DEPS)
      r_pb_dMdeps_p_even = r_pb_dMdeps_p_odd = (double)0.0;
      /* dMdeps_x = dM/deps0 M^{-1} gr */
      dMdeps_x = create_v_field();
      dslash_fn_field( M_inv_gr, dMdeps_x, EVENANDODD, fn_deps );
#endif

#ifdef CHEM_POT
      pb_dMdmu_p_e = pb_dMdmu_p_o = dcmplx((double)0.0,(double)0.0);
      pb_d2Mdmu2_p_e = pb_d2Mdmu2_p_o = dcmplx((double)0.0,(double)0.0);

      /* dM_M_inv_gr = dM/dmu * M_inv_gr */
      /* d2M_M_inv_gr = d2M/dmu2 * M_inv_gr */
      dM_M_inv_gr = create_v_field();
      chem_pot_tshift(fn, dM_M_inv_gr, M_inv_gr, 3., -1.);
      d2M_M_inv_gr = create_v_field();
      chem_pot_tshift(fn, d2M_M_inv_gr, M_inv_gr, 9., 1.);

#endif

      /* fermion action = M_gr.M_inv_gr */
      /* psi-bar-psi on even sites = gr.M_inv_gr */
      FOREVENFIELDSITES(i){
	rfaction += su3_rdot( M_gr+i, M_inv_gr+i );
	cc = su3_dot( gr+i, M_inv_gr+i );
	CSUM(pbp_e, cc);

#ifdef DM_DU0
	/* r_pb_dMdu_p_even = gr * dM/du0 M^{-1} gr |even*/
	r_pb_dMdu_p_even += su3_rdot( gr+i, dMdu_x+i );
#endif

#if FERM_ACTION == HISQ & defined(DM_DEPS)
	/* r_pb_dMdu_p_even = gr * dM/du0 M^{-1} gr |even*/
	r_pb_dMdeps_p_even += su3_rdot( gr+i, dMdeps_x+i );
#endif

#ifdef CHEM_POT
	/* Compute pb_dMdmu_p, pb_d2Mdmu2_p and dM_M_inv on even sites */
	cc = su3_dot( gr+i, dM_M_inv_gr+i);
	CSUM(pb_dMdmu_p_e, cc);
	cc = su3_dot( gr+i, d2M_M_inv_gr+i);
	CSUM(pb_d2Mdmu2_p_e, cc);
#endif
      }

      /* psi-bar-psi on odd sites */
      FORODDFIELDSITES(i){
	cc = su3_dot( gr+i, M_inv_gr+i );
	CSUM(pbp_o, cc);
#ifdef DM_DU0
	/* r_pb_dMdu_p_odd = gr * dM/du0 M^{-1} gr |odd*/
	r_pb_dMdu_p_odd += su3_rdot( gr+i, dMdu_x+i );
#endif

#if FERM_ACTION == HISQ & defined(DM_DEPS)
	/* r_pb_dMdu_p_odd = gr * dM/du0 M^{-1} gr |odd*/
	r_pb_dMdeps_p_odd += su3_rdot( gr+i, dMdeps_x+i );
#endif

#ifdef CHEM_POT
	/* Compute pb_dMdmu_P, pb_d2Mdmu2_p and dM_M_inv on odd sites */
	cc = su3_dot( gr+i, dM_M_inv_gr+i);
	CSUM(pb_dMdmu_p_o, cc);
	cc = su3_dot( gr+i, d2M_M_inv_gr+i);
	CSUM(pb_d2Mdmu2_p_o, cc);
#endif
      }

#ifdef CHEM_POT
      destroy_v_field(d2M_M_inv_gr); d2M_M_inv_gr = NULL;
#endif
      destroy_v_field(M_gr); M_gr = NULL;

      g_dcomplexsum( &pbp_o );
      g_dcomplexsum( &pbp_e );
      g_doublesum( &rfaction );
      
#ifdef DM_DU0
      destroy_v_field( dMdu_x ); dMdu_x = NULL;
      g_doublesum( &r_pb_dMdu_p_even );
      g_doublesum( &r_pb_dMdu_p_odd );
      r_pb_dMdu_p_even *= (2.0/(double)volume);
      r_pb_dMdu_p_odd *= (2.0/(double)volume);
      node0_printf("PB_DMDU_P: mass %e  %e  %e ( %d of %d )\n", mass,
		   r_pb_dMdu_p_even, r_pb_dMdu_p_odd, jpbp_reps+1, npbp_reps);
#endif

#if FERM_ACTION == HISQ & defined(DM_DEPS)
      destroy_v_field( dMdeps_x ); dMdeps_x = NULL;
      g_doublesum( &r_pb_dMdeps_p_even );
      g_doublesum( &r_pb_dMdeps_p_odd );
      r_pb_dMdeps_p_even *= (2.0/(double)volume);
      r_pb_dMdeps_p_odd *= (2.0/(double)volume);
      node0_printf("PB_DMDEPS_P: mass %e  %e  %e ( %d of %d )\n", mass,
		   r_pb_dMdeps_p_even, r_pb_dMdeps_p_odd, jpbp_reps+1, npbp_reps);
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
      /* Print results for pb_dMdmu_p and pb_d2Mdmu2_p */
      chem_pot_print1(pb_dMdmu_p_e, pb_dMdmu_p_o, pb_d2Mdmu2_p_e, pb_d2Mdmu2_p_o,
		      mass, jpbp_reps, npbp_reps);
#endif

#ifdef TR_MM_INV
      if(npbp_reps > 1){
	su3_vector *MM_inv_gr = create_v_field();
	double pbp_pbp = 0.0;

	mat_invert_uml_field( M_inv_gr, MM_inv_gr, qic, mass, fn );
	FORALLFIELDSITES(i){
	  pbp_pbp += su3_rdot( gr+i, MM_inv_gr+i );
	}
	g_doublesum( &pbp_pbp );
	pbp_pbp =  pbp_pbp*(1.0/(double)volume) ;
	node0_printf("TR_MM_INV: mass %e,  %e ( %d of %d )\n", mass,
		     pbp_pbp, jpbp_reps+1, npbp_reps);
	destroy_v_field(MM_inv_gr);
      }
#endif
      destroy_v_field(M_inv_gr); M_inv_gr = NULL;

#ifdef CHEM_POT
      /* M_inv_dM_M_inv_gr = M^{-1} dM_M_inv_gr */

      M_inv_dM_M_inv_gr = create_v_field();
      mat_invert_uml_field( dM_M_inv_gr, M_inv_dM_M_inv_gr, qic, mass, fn );
      destroy_v_field(dM_M_inv_gr); dM_M_inv_gr = NULL;

      /* dM_M_inv_dM_M_inv_gr = dM/dmu M_inv_dM_M_inv_gr */

      dM_M_inv_dM_M_inv_gr = create_v_field();
      chem_pot_tshift(fn, dM_M_inv_dM_M_inv_gr, M_inv_dM_M_inv_gr, 3., -1.);
      destroy_v_field(M_inv_dM_M_inv_gr); M_inv_dM_M_inv_gr = NULL;

      /* Compute MidM_MidM */
      MidM_MidM = (double)0.0;
      FORALLFIELDSITES(i){
	MidM_MidM += su3_rdot( gr+i, dM_M_inv_dM_M_inv_gr+i);
      }

      destroy_v_field(dM_M_inv_dM_M_inv_gr); dM_M_inv_dM_M_inv_gr = NULL;

      chem_pot_print2(MidM_MidM, jpbp_reps, npbp_reps, mass);

#endif
      destroy_v_field(gr); gr = NULL;

    } /* jpbp_reps */
}

/* Wrapper for obsolete call */

void f_meas_imp( int npbp_reps, int prec, 
		 field_offset phi_off, field_offset xxx_off, Real mass,
		 int naik_term_epsilon_index, fermion_links_t *fl){
  quark_invert_control qic;
  qic.prec       = prec;
  qic.min        = 0;
  qic.max        = niter;
  qic.nrestart   = nrestart;
  qic.parity     = EVENANDODD;
  qic.start_flag = 0;
  qic.nsrc       = 1;
  qic.resid      = sqrt(rsqprop);
  qic.relresid   = 0;

  f_meas_imp_field( npbp_reps, &qic, mass, naik_term_epsilon_index, fl );
}

/* Entry point for multiple masses.  Saves a few cycles because one
   inversion can be done with the multimass inverter */

void f_meas_imp_multi( int n_masses, int npbp_reps, quark_invert_control *qic, 
		       ks_param *ksp, fermion_links_t *fl){

  Real *mass = create_real_array(n_masses);
  imp_ferm_links_t **fn = get_fm_links(fl);

#ifdef DM_DU0
  imp_ferm_links_t *fn_du0 = get_fm_du0_links(fl)[0];
  su3_vector **dMdu_x = NULL;
#endif

#if FERM_ACTION == HISQ & defined(DM_DEPS)
  imp_ferm_links_t *fn_deps = get_fn_deps_links(fl);
  su3_vector **dMdeps_x = NULL;
#endif

  int i, j;
  int jpbp_reps;
  su3_vector *gr = NULL;
  su3_vector **M_gr = NULL;
  su3_vector **M_inv_gr = NULL;
  imp_ferm_links_t **fn_multi;

  /* Load masses from ks_param */
  for(j = 0; j < n_masses; j++)
    mass[j] = ksp[j].mass;

  /* Load pointers for fermion links, based on Naik epsilon indices */
  fn_multi = (imp_ferm_links_t **)malloc(sizeof(imp_ferm_links_t *)*n_masses);
  for(j = 0; j < n_masses; j++)
    fn_multi[j] = fn[ksp[j].naik_term_epsilon_index];

  for(jpbp_reps = 0; jpbp_reps < npbp_reps; jpbp_reps++){
      
    /* Make random source, and do inversion */
    /* generate gr random; M_gr = M gr */
    gr = create_v_field();
#ifndef Z2RSOURCE
    grsource_plain_field( gr, EVENANDODD );
#else
    z2rsource_plain_field( gr, EVENANDODD );
#endif

    M_gr = create_su3_vector_array(n_masses);
    for(j = 0; j < n_masses; j++)
      ks_dirac_adj_op( gr, M_gr[j], mass[j], EVENANDODD, fn_multi[j] );

    /* M_inv_gr = M^{-1} gr */

    M_inv_gr = create_su3_vector_array(n_masses);
    total_iters += mat_invert_multi( gr, M_inv_gr, ksp, n_masses, qic, fn_multi );
      
#ifdef DM_DU0
    /* dMdu_x = dM/du0 M^{-1} gr */
    dMdu_x = create_su3_vector_array(n_masses);
    for(j = 0; j < n_masses; j++){
      dslash_fn_field( M_inv_gr[j], dMdu_x[j], EVENANDODD, fn_du0 );
    }
#endif

#if FERM_ACTION == HISQ & defined(DM_DEPS)
    /* dMdeps_x = dM/deps M^{-1} gr */
    dMdeps_x = create_su3_vector_array(n_masses);
    for(j = 0; j < n_masses; j++){
      dslash_fn_field( M_inv_gr[j], dMdeps_x[j], EVENANDODD, fn_deps );
    }
#endif

#ifdef CHEM_POT
#endif

    for(j = 0; j < n_masses; j++){
      complex cc;
      double rfaction = 0.0;
      double_complex pbp_e = dcmplx(0.0,0.0);
      double_complex pbp_o = dcmplx(0.0,0.0);
      Real r_psi_bar_psi_even, i_psi_bar_psi_even;
      Real r_psi_bar_psi_odd, i_psi_bar_psi_odd;
      Real r_ferm_action;

#ifdef DM_DU0
      double r_pb_dMdu_p_even = 0.0;
      double r_pb_dMdu_p_odd = 0.0;
#endif

#if FERM_ACTION == HISQ & defined(DM_DEPS)
      double r_pb_dMdeps_p_even = 0.0;
      double r_pb_dMdeps_p_odd = 0.0;
#endif

#ifdef CHEM_POT
      su3_vector *dM_M_inv_gr = NULL;
      su3_vector *d2M_M_inv_gr = NULL;
      su3_vector *M_inv_dM_M_inv_gr = NULL;
      su3_vector *dM_M_inv_dM_M_inv_gr = NULL;
      double_complex pb_dMdmu_p_e = dcmplx(0.0,0.0);
      double_complex pb_dMdmu_p_o = dcmplx(0.0,0.0);
      double_complex pb_d2Mdmu2_p_e = dcmplx(0.0,0.0);
      double_complex pb_d2Mdmu2_p_o = dcmplx(0.0,0.0);
      double MidM_MidM = 0.0;

      /* dM_M_inv_gr = dM/dmu * M_inv_gr */
      /* d2M_M_inv_gr = d2M/dmu2 * M_inv_gr */
      dM_M_inv_gr = create_v_field();
      d2M_M_inv_gr = create_v_field();
      chem_pot_tshift(fn_multi[j], dM_M_inv_gr, M_inv_gr[j], 3., -1.);
      chem_pot_tshift(fn_multi[j], d2M_M_inv_gr, M_inv_gr[j], 9., 1.);

#endif

      /* fermion action = M_gr.M_inv_gr */
      /* psi-bar-psi on even sites = gr.M_inv_gr */
      FOREVENFIELDSITES(i){
	rfaction += su3_rdot( M_gr[j]+i, M_inv_gr[j]+i );
	cc = su3_dot( gr+i, M_inv_gr[j]+i );
	CSUM(pbp_e, cc);

#ifdef DM_DU0
	/* r_pb_dMdu_p_even = gr * dM/du0 M^{-1} gr |even*/
	r_pb_dMdu_p_even += su3_rdot( gr+i, dMdu_x[j]+i );
#endif

#if FERM_ACTION == HISQ & defined(DM_DEPS)
	/* r_pb_dMdeps_p_even = gr * dM/deps0 M^{-1} gr |even*/
	r_pb_dMdeps_p_even += su3_rdot( gr+i, dMdeps_x[j]+i );
#endif

#ifdef CHEM_POT
	/* Compute pb_dMdmu_p, pb_d2Mdmu2_p and dM_M_inv on even sites */
	cc = su3_dot( gr+i, dM_M_inv_gr+i);
	CSUM(pb_dMdmu_p_e, cc);
	cc = su3_dot( gr+i, d2M_M_inv_gr+i);
	CSUM(pb_d2Mdmu2_p_e, cc);
#endif
      }

      /* psi-bar-psi on odd sites */
      FORODDFIELDSITES(i){
	cc = su3_dot( gr+i, M_inv_gr[j]+i );
	CSUM(pbp_o, cc);
#ifdef DM_DU0
	/* r_pb_dMdu_p_odd = gr * dM/du0 M^{-1} gr |odd*/
	r_pb_dMdu_p_odd += su3_rdot( gr+i, dMdu_x[j]+i );
#endif

#if FERM_ACTION == HISQ & defined(DM_DEPS)
	/* r_pb_dMdeps_p_odd = gr * dM/deps0 M^{-1} gr |odd*/
	r_pb_dMdeps_p_odd += su3_rdot( gr+i, dMdeps_x[j]+i );
#endif

#ifdef CHEM_POT
	/* Compute pb_dMdmu_P, pb_d2Mdmu2_p and dM_M_inv on odd sites */
	cc = su3_dot( gr+i, dM_M_inv_gr+i);
	CSUM(pb_dMdmu_p_o, cc);
	cc = su3_dot( gr+i, d2M_M_inv_gr+i);
	CSUM(pb_d2Mdmu2_p_o, cc);
#endif
      }

#ifdef CHEM_POT
      destroy_v_field(d2M_M_inv_gr); d2M_M_inv_gr = NULL;
#endif
      destroy_v_field(M_gr[j]); M_gr[j] = NULL;

      g_dcomplexsum( &pbp_o );
      g_dcomplexsum( &pbp_e );
      g_doublesum( &rfaction );
      
#ifdef DM_DU0
      destroy_v_field(dMdu_x[j]); dMdu_x[j] = NULL;
      g_doublesum( &r_pb_dMdu_p_even );
      g_doublesum( &r_pb_dMdu_p_odd );
      r_pb_dMdu_p_even *= (2.0/(double)volume);
      r_pb_dMdu_p_odd *= (2.0/(double)volume);
      node0_printf("PB_DMDU_P: mass %e  %e  %e ( %d of %d )\n", mass[j],
		   r_pb_dMdu_p_even, r_pb_dMdu_p_odd, jpbp_reps+1, npbp_reps);
#endif

#if FERM_ACTION == HISQ & defined(DM_DEPS)
      destroy_v_field(dMdeps_x[j]); dMdeps_x[j] = NULL;
      g_doublesum( &r_pb_dMdeps_p_even );
      g_doublesum( &r_pb_dMdeps_p_odd );
      r_pb_dMdeps_p_even *= (2.0/(double)volume);
      r_pb_dMdeps_p_odd *= (2.0/(double)volume);
      node0_printf("PB_DMDEPS_P: mass %e  %e  %e ( %d of %d )\n", mass[j],
		   r_pb_dMdeps_p_even, r_pb_dMdeps_p_odd, jpbp_reps+1, npbp_reps);
#endif

      r_psi_bar_psi_odd =  pbp_o.real*(2.0/(double)volume) ;
      i_psi_bar_psi_odd =  pbp_o.imag*(2.0/(double)volume) ;
      r_psi_bar_psi_even =  pbp_e.real*(2.0/(double)volume) ;
      i_psi_bar_psi_even =  pbp_e.imag*(2.0/(double)volume) ;
      r_ferm_action =  rfaction*(1.0/(double)volume) ;
      node0_printf("PBP: mass %e     %e  %e  %e  %e ( %d of %d )\n", mass[j],
		   r_psi_bar_psi_even, r_psi_bar_psi_odd,
		   i_psi_bar_psi_even, i_psi_bar_psi_odd,
		   jpbp_reps+1, npbp_reps);
      node0_printf("FACTION: mass = %e,  %e ( %d of %d )\n", mass[j],
		   r_ferm_action, jpbp_reps+1, npbp_reps);
      
#ifdef CHEM_POT
      /* Print results for pb_dMdmu_p and pb_d2Mdmu2_p */
      chem_pot_print1(pb_dMdmu_p_e, pb_dMdmu_p_o, pb_d2Mdmu2_p_e, pb_d2Mdmu2_p_o,
		      mass[j], jpbp_reps, npbp_reps);
#endif
      
#ifdef TR_MM_INV
      if(npbp_reps > 1){
	su3_vector *MM_inv_gr = create_v_field();
	double pbp_pbp = 0.0;
	
	mat_invert_uml_field( M_inv_gr[j], MM_inv_gr, qic, mass[j], fn_multi[j] );
	FORALLFIELDSITES(i){
	  pbp_pbp += su3_rdot( gr+i, MM_inv_gr+i );
	}
	g_doublesum( &pbp_pbp );
	pbp_pbp =  pbp_pbp*(1.0/(double)volume) ;
	node0_printf("TR_MM_INV: mass %e,  %e ( %d of %d )\n", mass[j],
		     pbp_pbp, jpbp_reps+1, npbp_reps);
	destroy_v_field(MM_inv_gr); MM_inv_gr = NULL;
      }
#endif
      destroy_v_field(M_inv_gr[j]); M_inv_gr[j] = NULL;
      
#ifdef CHEM_POT
      /* M_inv_dM_M_inv_gr = M^{-1} dM_M_inv_gr */
      
      M_inv_dM_M_inv_gr = create_v_field();
      mat_invert_uml_field( dM_M_inv_gr, M_inv_dM_M_inv_gr, qic, mass[j], fn_multi[j] );
      destroy_v_field(dM_M_inv_gr); dM_M_inv_gr = NULL;
      
      /* dM_M_inv_dM_M_inv_gr = dM/dmu M_inv_dM_M_inv_gr */
      
      dM_M_inv_dM_M_inv_gr = create_v_field();
      chem_pot_tshift(fn_multi[j], dM_M_inv_dM_M_inv_gr, M_inv_dM_M_inv_gr, 3., -1.);
      destroy_v_field(M_inv_dM_M_inv_gr); M_inv_dM_M_inv_gr = NULL;
      
      /* Compute MidM_MidM */
      FORALLFIELDSITES(i){
	MidM_MidM += su3_rdot( gr+i, dM_M_inv_dM_M_inv_gr+i);
      }
      
      destroy_v_field(dM_M_inv_dM_M_inv_gr); dM_M_inv_dM_M_inv_gr = NULL;
      
      chem_pot_print2(MidM_MidM, jpbp_reps, npbp_reps, mass[j]);
      
#endif
      
    } /* for j masses */
    
    destroy_v_field(gr); gr = NULL;
#ifdef DM_DU0
    destroy_su3_vector_array(dMdu_x,n_masses);
#endif
#if FERM_ACTION == HISQ & defined(DM_DEPS)
    destroy_su3_vector_array(dMdeps_x,n_masses);
#endif
    destroy_su3_vector_array(M_inv_gr,n_masses);
    destroy_su3_vector_array(M_gr,n_masses);
    
  } /* jpbp_reps */
		     
  free(fn_multi);
  destroy_real_array(mass); mass = NULL;
}

