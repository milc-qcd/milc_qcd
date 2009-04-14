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

void update_h_rhmc( Real eps, su3_vector **multi_x ){
#ifdef FN
  free_fn_links(&fn_links);
  free_fn_links(&fn_links_dmdu0);
#endif
  /*  node0_printf("update_h_rhmc:\n"); */
  /* gauge field force */
  rephase(OFF);
  imp_gauge_force(eps,F_OFFSET(mom));
  rephase(ON);
  /* fermionic force */
  
  update_h_fermion( eps,  multi_x );
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
void update_h_fermion( Real eps, su3_vector **multi_x ){
  int iphi,jphi;
  Real final_rsq;
  int i,j;
  int order, tmporder;
  Real *residues,*allresidues;
  Real *roots;

  //node0_printf("update_h_rhmc: EXPERIMENTAL force call\n");
  //fflush(stdout);

  /* Algorithm sketch: assemble multi_x with all |X> fields,
     then call force routine for each part (so far we have to parts:
     zero correction to Naik and non-zero correction to Naik */

  allresidues = (Real *)malloc(n_order_naik_total*sizeof(Real));

  // Group the fermion force calculation according to sets of like
  // path coefficients.
  tmporder = 0;
  iphi = 0;
  for( i=0; i<n_naiks; i++ ) {
    for( jphi=0; jphi<n_pseudo_naik[i]; jphi++ ) {
#ifdef HISQ
      fn_links.hl.current_X_set = i; // which X set we need
#endif
      load_ferm_links(&fn_links, &ks_act_paths);

      // Add the current pseudofermion to the current set
      order = rparam[iphi].MD.order;
      residues = rparam[iphi].MD.res;
      roots = rparam[iphi].MD.pole;

      // Compute ( M^\dagger M)^{-1} in xxx_even
      // Then compute M*xxx in temporary vector xxx_odd 
      /* See long comment at end of file */
	/* The diagonal term in M doesn't matter */
      ks_ratinv( F_OFFSET(phi[iphi]), multi_x+tmporder, roots, order, 
	         niter_md[iphi], rsqmin_md[iphi], prec_md[iphi], EVEN, 
	         &final_rsq, &fn_links );

      for(j=0;j<order;j++){
	dslash_field( multi_x[tmporder+j], multi_x[tmporder+j],  ODD,
		      &fn_links);
	allresidues[tmporder+j] = residues[j+1];
	// remember that residues[0] is constant, no force contribution.
      }
    tmporder += order;
    iphi++;
    }
  }

#ifdef MILC_GLOBAL_DEBUG
  node0_printf("update_h_rhmc: MULTI_X ASSEMBLED\n");fflush(stdout);
  node0_printf("update_h_rhmc: n_distinct_Naik=%d\n",n_naiks);
  for(j=0;j<n_naiks;j++)
    node0_printf("update_h_rhmc: orders[%d]=%d\n",j,n_orders_naik[j]);
  for(j=0;j<n_naiks;j++)
    node0_printf("update_h_rhmc: masses_Naik[%d]=%f\n",j,masses_naik[j]);
  fflush(stdout);
#endif /* MILC_GLOBAL_DEBUG */

  eo_fermion_force_multi( eps, allresidues, multi_x,
         n_order_naik_total, prec_ff, &fn_links, &ks_act_paths );

  free(allresidues);
} /* update_h_fermion */

