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
  int iphi;
  Real final_rsq;
  int j;
  int order, tmporder;
  int path_coeff_changed;
  Real *residues,*allresidues;
  Real *roots;
  
  for(tmporder=0,iphi = 0; iphi < n_pseudo; iphi++)
    tmporder+=rparam[iphi].MD.order;
  allresidues = (Real *)malloc(tmporder*sizeof(Real));

  load_ferm_links(&fn_links, &ks_act_paths);
  
  // Group the fermion force calculation according to sets of like
  // path coefficients.
  tmporder = 0;
  for(iphi = 0; iphi < n_pseudo; iphi++){
    // Remake the path tables if the coeffs change for this mass
    path_coeff_changed = make_path_table(&ks_act_paths, &ks_act_paths_dmdu0,
					 rparam[iphi].naik_term_mass);
    if(path_coeff_changed){
      // Invalidate only fat and long links and remake them
      invalidate_fn_links(&fn_links);
      load_ferm_links(&fn_links, &ks_act_paths);
      // Calculate the contribution of the previous set to the fermion force
      if(tmporder > 0){
	eo_fermion_force_multi( eps, allresidues, multi_x, 
				tmporder, prec_ff, &fn_links, 
				&ks_act_paths );
	tmporder = 0;
      }
    }
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
  }

  // Do the last set
  if(tmporder > 0)
    eo_fermion_force_multi( eps, allresidues, multi_x, tmporder, 
			    prec_ff, &fn_links, &ks_act_paths );

  free(allresidues);
} /* update_h_fermion */

