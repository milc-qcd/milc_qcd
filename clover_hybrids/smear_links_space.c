/************************** smear_links_space.c *******************************/
/* MIMD version 7 */
/* Doug Toussaint 1/30/96 */
/* 7/17/96 version4 to smear only in space directions */
/* 11/16/01 CD changed to call standardized ape_smear_dir */

/* Construct smeared links in dest[dir1].   Smeared link is sum
   of link plus all staples:

		---------		-------> dir1
		|	|
		|	|
		|	|		^
		X--------		| dir2
		|	|		|
		|	|
		|	|
		---------
   dest is probably secretly defined to overlap conjugate gradient
   vectors in the Wilson code.

   For the moment, we add the simple link in with relative weight one,
   arbitrary normalization factor to keep F_mu_nu reasonable size.
   */

/* SPACE_NORM_FACTOR = 1 / (SIMPLE_WEIGHT+3) is right for ordered
   lattice.  Something larger for disordered.  Perhaps
   1  / ( SIMPLE_WEIGHT + 3*u_0^2 ), where u_0 is tadpole
   improvement parameter, or perhaps determine NORM_FACTOR
   by trial at beginning of run. */

#include "cl_hyb_includes.h"

void smear_links( field_offset src, field_offset dest ){
  int dir1;
  Real simple_weight, norm_factor;
  Real staple_weight, link_u0;
  
  if(this_node==0){
    printf(
	   "Smearing: space_simple_weight = %e , space_norm_factor = %e\n",
	   space_simple_weight,space_norm_factor);
    printf(
	   "           time_simple_weight = %e ,  time_norm_factor = %e\n",
	   time_simple_weight,time_norm_factor);
  }
  
  for(dir1=XUP;dir1<=TUP;dir1++){
    if(dir1==TUP){
      simple_weight=time_simple_weight;
      norm_factor=time_norm_factor;
      link_u0 = 
	sqrt((1.0 - norm_factor*simple_weight)/(norm_factor*6.0));
    } else {
      simple_weight=space_simple_weight;
      norm_factor=space_norm_factor;
      link_u0 = 
	sqrt((1.0 - norm_factor*simple_weight)/(norm_factor*4.0));
    }
    staple_weight = 1/simple_weight;
    
    /* Smear links in the specified direction with specified
       weights using space only and no reprojection onto SU(3) */ 
    
    ape_smear_dir( src, dir1, dest+sizeof(su3_matrix)*dir1, 
		   staple_weight, link_u0, 1, 0, 0.);
    
  } /*dir1 loop */
} /* smear_links */

