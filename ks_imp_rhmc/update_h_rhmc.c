/****** update_h_rhmc.c  -- ******************/
/* MIMD version 7 */
/* updates momentum matrices for improved action  with RHMC algorith*/
/* D.T. & J.H., naik term    8/96
*  D.T., fat link fermion term 5/97
*  D.T. general quark action 1/98
*  D.T. two types of quarks 3/99
*  T.D. and A.H. improved gauge updating spliced in 5/97
*  D.T. first try at RHMC version 12/05
*  D.T. 3/07 Gang together multiple pseudofermion terms in update_h_fermion
*/

#include "ks_imp_includes.h"	/* definitions files and prototypes */

int update_h_rhmc( Real eps, su3_vector **multi_x ){
  int iters;
#ifdef FN
  invalidate_fermion_links(fn_links);
#endif
  /*  node0_printf("update_h_rhmc:\n"); */
  /* gauge field force */
  rephase(OFF);
#ifndef U1_ONLY
  imp_gauge_force_ks(eps,F_OFFSET(mom));
#endif /* U1_ONLY */
#ifdef HAVE_U1
  gauge_force_u1(eps);
#endif /* HAVE_U1 */
  rephase(ON);
#ifndef PURE_GAUGE_U1
  /* fermionic force */
  iters = update_h_fermion( eps,  multi_x );
  return iters;
#else
  return 0;
#endif /* PURE_GAUGE_U1 */
} /* update_h_rhmc */

// gauge and fermion force parts separately, for algorithms that use
// different time steps for them
void update_h_gauge( Real eps ){
  /* node0_printf("update_h_gauge:\n");*/
  /* gauge field force */
  imp_gauge_force_ks(eps,F_OFFSET(mom));
  rephase(OFF);
#ifndef U1_ONLY
  imp_gauge_force_ks(eps,F_OFFSET(mom));
#endif /* U1_ONLY */
#ifdef HAVE_U1
  gauge_force_u1(eps);
#endif
  rephase(ON);
} /* update_h_gauge */

// fermion force update grouping pseudofermions with the same path coeffs
int update_h_fermion( Real eps, su3_vector **multi_x ){
#ifndef PURE_GAUGE_U1
  int iphi,jphi;
  Real final_rsq;
  int i,j,n;
  int order, tmporder;
  Real *residues,*allresidues;
  Real *roots;
  int iters = 0;
  imp_ferm_links_t *fn;

  /* Algorithm sketch: assemble multi_x with all |X> fields,
     then call force routine for each part (so far we have to parts:
     zero correction to Naik and non-zero correction to Naik */

  allresidues = (Real *)malloc(n_order_naik_total*sizeof(Real));

  // Group the fermion force calculation according to sets of like
  // path coefficients.
  tmporder = 0;
  iphi = 0;
  restore_fermion_links_from_site(fn_links, prec_md[0]);
#if FERM_ACTION == HISQ
  n = fermion_links_get_n_naiks(fn_links);
//printf("update_h_rhmc fermion_links_get_n_naiks %d\n", n);
#else
  n = 1;
#endif
  for( i=0; i<n; i++ ) {
    /* Assume prec_md is the same for all pseudo-fermions
       with the same Naik epsilon */
    fn = get_fm_links(fn_links, i);
    for( jphi=0; jphi<n_pseudo_naik[i]; jphi++ ) {
      
#ifdef HAVE_U1
      /* can be improved by checking whether the charge is changed */
      current_charge_u1 = 1.0 * pseudo_charges[iphi];
      u1phase_on(current_charge_u1, u1_A);
      invalidate_fermion_links(fn_links);
#endif
      restore_fermion_links_from_site(fn_links, prec_md[iphi]);
      fn = get_fm_links(fn_links, i);

      // Add the current pseudofermion to the current set
      order = rparam[iphi].MD.order;
      residues = rparam[iphi].MD.res;
      roots = rparam[iphi].MD.pole;

      // Compute ( M^\dagger M)^{-1} in xxx_even
      // Then compute M*xxx in temporary vector xxx_odd 
      /* See long comment at end of file */
	/* The diagonal term in M doesn't matter */
      iters += ks_ratinv( F_OFFSET(phi[iphi]), multi_x+tmporder, roots, residues,
                          order, niter_md[iphi], rsqmin_md[iphi], prec_md[iphi], EVEN,
			  &final_rsq, fn, i, rparam[iphi].naik_term_epsilon );

      for(j=0;j<order;j++){
	dslash_field( multi_x[tmporder+j], multi_x[tmporder+j],  ODD, fn);
	allresidues[tmporder+j] = residues[j+1];
	// remember that residues[0] is constant, no force contribution.
      }
#ifdef HAVE_U1
      /* Unapply the U(1) field phases */
      u1phase_off();
      invalidate_fermion_links(fn_links);
#endif
    tmporder += order;
    iphi++;
    }
    destroy_fn_links(fn);
  }

#ifdef MILC_GLOBAL_DEBUG
  node0_printf("update_h_rhmc: MULTI_X ASSEMBLED\n");fflush(stdout);
  node0_printf("update_h_rhmc: n_distinct_Naik=%d\n",n);
  for(j=0;j<n;j++)
    node0_printf("update_h_rhmc: orders[%d]=%d\n",j,n_orders_naik[j]);
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  for(j=0;j<n;j++)
//    node0_printf("update_h_rhmc: masses_Naik[%d]=%f\n",j,fn_links.hl.eps_naik[j]);
#endif
  fflush(stdout);
#endif /* MILC_GLOBAL_DEBUG */

#ifdef HAVE_U1
  tmporder = 0;
  /* loop over different unique charges */
  for ( i = 0; i < n_charges_uniq; i++ ) {
    /* current_charge_u1 = 1.0*charges_uniq[i]; */
#ifdef U1_DEBUG
    node0_printf("update_h_rhmc.c i, charges_uniq[i]: %d %e\n", i,
                 charges_uniq[i]);
    node0_printf("update_h_rhmc.c i, tmporder, n_orders_naik_charge[i], "
                 "n_orders_naik_charge_heavy[i]: %d %d %d %d\n",
                 i, tmporder, n_orders_naik_charge[i],
                 n_orders_naik_charge_heavy[i]);
#endif
    /* u1phase_on(current_charge_u1, u1_A); */
    u1phase_on(charges_uniq[i], u1_A);
    invalidate_fermion_links(fn_links);
    restore_fermion_links_from_site(fn_links, prec_ff);
    eo_fermion_force_multi_su3_u1(eps, allresidues + tmporder,
                                  multi_x + tmporder,
                                  n_orders_naik_charge[i],
                                  n_orders_naik_charge_heavy[i],
                                  prec_ff, fn_links, charges_uniq[i]);
#ifdef U1_DEBUG
    node0_printf("update_h_rhmc.c after i, charges_uniq[i]: %d %e\n", i,
                 charges_uniq[i]);
#endif
    u1phase_off();
    invalidate_fermion_links(fn_links);
    tmporder += n_orders_naik_charge[i];
  }
#else
  restore_fermion_links_from_site(fn_links, prec_ff);
  eo_fermion_force_multi( eps, allresidues, multi_x,
                          n_order_naik_total, prec_ff, fn_links );
#endif

  free(allresidues);
  return iters;
#else
  return 0;
#endif /* PURE_GAUGE_U1 */
} /* update_h_fermion */

