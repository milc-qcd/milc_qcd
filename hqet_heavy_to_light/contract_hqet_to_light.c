/************* contract_hqet_to_light.c **************************/
/* MIMD version 4 */
/* 
 *  Do the Wick contractions for the hqet to light 
 *  three point function.
 *
 */

/* Modifications

   C. DeTar 4/30/97 Now working from a table of source-sink gamma matrix pairs
   C. DeTar 5/24/97 Now computes form factor with SW rotation

 */


#include "hqet_light_includes.h"
#include "opertypes.h"
#include "corrlist.h"

#ifdef DEBUGDEF
#include DEBUGDEF
#endif


/*  Do the Wick contractions of the quark propagators to form the 
 *  three point functions
 *
 *
 *   Function arguments
 * 
 *    On input
 *        q_zonked      :: site structure pointer to spin_wilson_vector 
 *                         [[ zonked quark ]]
 *        q_zonked_rot  :: site structure pointer to spin_wilson_vector 
 *                         [[ SW rotated zonked quark ]]
 *        q_sequential  :: site structure pouinter to spin_wilson_vector 
 *                         [[ sequential inverted quark ]]
 *        zonked_pt :: which zonked quark 
 *        spect_pt  :: which spectator quark
 *
 *    On output
 *        corr :: the three point function correlators 
 *
 *  
 *   This code calls the general purpose function  "meson_cont_mom",
 *   which does the low level contraction. This routine includes the 
 *   hermitian conjugate and multiplies the zonked quark by gamma_5
 *   on both dides.
 */



void contract_hqet_to_light(complex *corr, 
			    field_offset q_zonked, 
			    field_offset q_zonked_rot,
			    field_offset q_sequential,
			    int vel_pt,  int zonked_pt, int spect_pt)
{
  double t_start ;
  int base_pt, q_stride, op_stride ;

  t_start = dclock() ;

  /* Compute partial offset for storage of result in corr[] */

  base_pt    = HQET_FORM_WHERE(0,zonked_pt,spect_pt,0,vel_pt, 0 ) ; 
  q_stride   = HQET_FORM_WHERE(0,0,        0,       1,0,      0 ) ;
  op_stride  = HQET_FORM_WHERE(0,0,        0,       0,0,      1 ) ;

  /* First, contract zonked and sequential */

  meson_cont_mom(corr , q_zonked, q_sequential, 
		 base_pt, q_stride, op_stride, 
		 hqet_to_light, MAX_THREEPT) ;

  /* Second, contract rotated zonked and sequential 
     Results go to second half of corr */

  base_pt += op_stride*MAX_THREEPT;
  meson_cont_mom(corr , q_zonked_rot, q_sequential, 
		 base_pt, q_stride, op_stride, 
		 hqet_to_light, MAX_THREEPT) ;

  IF_VERBOSE_ON(1)
    printf("contract_hqet_to_light::Time to Wick contract hqet-->light correlators = %g sec\n",dclock() - t_start) ;
  
} 
