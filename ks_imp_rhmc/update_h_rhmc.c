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
  free_fn_links();
#endif
  /*  node0_printf("update_h_rhmc:\n"); */
  /* gauge field force */
  rephase(OFF);
  imp_gauge_force(eps,F_OFFSET(mom));
  rephase(ON);
  /* fermionic force */
  
  update_h_fermion( eps,  multi_x );
  //for(iphi = 0; iphi < n_pseudo; iphi++){
    //eo_fermion_force_rhmc( eps, &rparam[iphi].MD, 
			   //multi_x, F_OFFSET(phi[iphi]), 
			   //rsqmin_md[iphi], niter_md[iphi], prec_md[iphi] );
  //}
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

// void update_h_fermion( Real eps, su3_vector **multi_x ){
//   int iphi;
// #ifdef FN
//   free_fn_links();
// #endif
//   /*node0_printf("update_h_fermion:\n"); */
//   /* fermionic force */
//   
//   for(iphi = 0; iphi < n_pseudo; iphi++){
//     eo_fermion_force_rhmc( eps, &rparam[iphi].MD, 
// 			   multi_x, F_OFFSET(phi[iphi]),
// 			   rsqmin_md[iphi], niter_md[iphi], prec_md[iphi]);
//   }
// } /* update_h_fermion */

// fermion force update combining all pseudofermions into single force call
void update_h_fermion( Real eps, su3_vector **multi_x ){
  int iphi;
#ifdef FN
  free_fn_links();
#endif
  /*node0_printf("update_h_fermion:\n"); */
  /* fermionic force */
    Real final_rsq;
    int j;
    int order,totalorder,tmporder;
    Real *residues,*allresidues;
    Real *roots;
  
  for(totalorder=0,iphi = 0; iphi < n_pseudo; iphi++)totalorder+=rparam[iphi].MD.order;
  allresidues = (Real *)malloc(totalorder*sizeof(Real));

  for(tmporder=0,iphi = 0; iphi < n_pseudo; iphi++){
    order = rparam[iphi].MD.order;
    residues = rparam[iphi].MD.res;
    roots = rparam[iphi].MD.pole;

    // Compute ( M^\dagger M)^{-1} in xxx_even
    // Then compute M*xxx in temporary vector xxx_odd 
    /* See long comment at end of file */
	/* The diagonal term in M doesn't matter */
    ks_ratinv( F_OFFSET(phi[iphi]), &(multi_x[tmporder]), roots, order, niter_md[iphi],
	       rsqmin_md[iphi], prec_md[iphi], EVEN, &final_rsq );

    for(j=0;j<order;j++){
	dslash_field( multi_x[tmporder+j], multi_x[tmporder+j],  ODD );
	allresidues[tmporder+j] = residues[j+1];
	// remember that residues[0] is constant, no force contribution.
    }

    tmporder += order;
  }
  eo_fermion_force_multi( eps, allresidues, multi_x, totalorder, prec_ff );
  free(allresidues);
} /* update_h_fermion */

