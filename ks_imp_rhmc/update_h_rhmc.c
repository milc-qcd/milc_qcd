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
  imp_gauge_force(eps,F_OFFSET(mom));
  rephase(ON);
  /* fermionic force */
  
  iters = update_h_fermion( eps,  multi_x );
  return iters;
} /* update_h_rhmc */

// gauge and fermion force parts separately, for algorithms that use
// different time steps for them
void update_h_gauge( Real eps ){
  /* node0_printf("update_h_gauge:\n");*/
  /* gauge field force */
  rephase(OFF);
  imp_gauge_force(eps,F_OFFSET(mom));
  rephase(ON);
} /* update_h_gauge */

// fermion force update grouping pseudofermions with the same path coeffs
int update_h_fermion( Real eps, su3_vector **multi_x ){
  int iphi,jphi;
  Real final_rsq;
  int i,j,n;
  int order, tmporder;
  Real *residues,*allresidues;
  Real *roots;
  int iters = 0;
  imp_ferm_links_t **fn;

  /* Algorithm sketch: assemble multi_x with all |X> fields,
     then call force routine for each part (so far we have to parts:
     zero correction to Naik and non-zero correction to Naik */

  allresidues = (Real *)malloc(n_order_naik_total*sizeof(Real));

  // Group the fermion force calculation according to sets of like
  // path coefficients.
  tmporder = 0;
  iphi = 0;
#if FERM_ACTION == HISQ
  n = fermion_links_get_n_naiks(fn_links);
#else
  n = 1;
#endif
  for( i=0; i<n; i++ ) {
    for( jphi=0; jphi<n_pseudo_naik[i]; jphi++ ) {
      restore_fermion_links_from_site(fn_links, prec_md[iphi]);
      fn = get_fm_links(fn_links);

      // Add the current pseudofermion to the current set
      order = rparam[iphi].MD.order;
      residues = rparam[iphi].MD.res;
      roots = rparam[iphi].MD.pole;

      // Compute ( M^\dagger M)^{-1} in xxx_even
      // Then compute M*xxx in temporary vector xxx_odd 
      /* See long comment at end of file */
	/* The diagonal term in M doesn't matter */
      iters += ks_ratinv( F_OFFSET(phi[iphi]), multi_x+tmporder, roots, order, 
			  niter_md[iphi], rsqmin_md[iphi], prec_md[iphi], EVEN, 
			  &final_rsq, fn[i], 
			  i, rparam[iphi].naik_term_epsilon );

      for(j=0;j<order;j++){
	dslash_field( multi_x[tmporder+j], multi_x[tmporder+j],  ODD,
		      fn[i]);
	allresidues[tmporder+j] = residues[j+1];
	// remember that residues[0] is constant, no force contribution.
      }
    tmporder += order;
    iphi++;
    }
  }

#ifdef MILC_GLOBAL_DEBUG
  node0_printf("update_h_rhmc: MULTI_X ASSEMBLED\n");fflush(stdout);
  node0_printf("update_h_rhmc: n_distinct_Naik=%d\n",n);
  for(j=0;j<n;j++)
    node0_printf("update_h_rhmc: orders[%d]=%d\n",j,n_orders_naik[j]);
#if FERM_ACTION == HISQ
  for(j=0;j<n;j++)
    node0_printf("update_h_rhmc: masses_Naik[%d]=%f\n",j,fn_links.hl.eps_naik[j]);
#endif
  fflush(stdout);
#endif /* MILC_GLOBAL_DEBUG */

  restore_fermion_links_from_site(fn_links, prec_ff);
  eo_fermion_force_multi( eps, allresidues, multi_x,
			  n_order_naik_total, prec_ff, fn_links );

  free(allresidues);
  return iters;
} /* update_h_fermion */

