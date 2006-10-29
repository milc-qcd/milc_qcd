/****** update_h_rhmc.c  -- ******************/
/* updates momentum matrices for improved action  with RHMC algorith*/
/* D.T. & J.H., naik term    8/96
*  D.T., fat link fermion term 5/97
*  D.T. general quark action 1/98
*  D.T. two types of quarks 3/99
*  T.D. and A.H. improved gauge updating spliced in 5/97
*  D.T. first try at RHMC version 12/05
*  alg_flag needed for algorithms with different forces at different steps
*
* MIMD version 7 */

#include "ks_imp_includes.h"	/* definitions files and prototypes */

void update_h_rhmc( int alg_flag, Real eps, su3_vector **multi_x ){
  int iphi;
#ifdef FN
  free_fn_links();
#endif
  node0_printf("update_h_rhmc: alg_flag=%d\n",alg_flag);
  /* gauge field force */
  rephase(OFF);
  imp_gauge_force(eps,F_OFFSET(mom));
  rephase(ON);
  /* fermionic force */
  
  for(iphi = 0; iphi < nphi; iphi++){
    eo_fermion_force_rhmc( alg_flag, eps, &rparam[iphi].MD, 
			   multi_x, F_OFFSET(phi[iphi]) );
  }
} /* update_h_rhmc */

// gauge and fermion force parts separately, for algorithms that use
// different time steps for them
void update_h_gauge( int alg_flag, Real eps ){
  node0_printf("update_h_gauge: alg_flag=%d\n",alg_flag);
  /* gauge field force */
  rephase(OFF);
  imp_gauge_force(eps,F_OFFSET(mom));
  rephase(ON);
} /* update_h_gauge */

void update_h_fermion( int alg_flag, Real eps, su3_vector **multi_x ){
  int iphi;
#ifdef FN
  free_fn_links();
#endif
  node0_printf("update_h_fermion: alg_flag=%d\n",alg_flag);
  /* fermionic force */
  
  for(iphi = 0; iphi < nphi; iphi++){
    eo_fermion_force_rhmc( alg_flag, eps, &(rparam[iphi].MD), 
			   multi_x, F_OFFSET(phi[iphi]) );
  }
} /* update_h_fermion */
