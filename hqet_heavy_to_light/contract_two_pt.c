/************* contract_two_pt.c **************************/
/* MIMD version 4  ***/
/*  
 *  Do the Wick contractions for the relevant two point functions
 *
 *  This routine will need to be updated, when the light quark
 *  smearing is implemented properly
 */

/* Modifications

   C. McNeile 1997   Original
   C. DeTar 5/24/97  Correlator table now in corrlist.h

 */

#include "hqet_light_includes.h"
#include "opertypes.h"
#include "corrlist.h"

#ifdef DEBUGDEF
#include DEBUGDEF
#endif

/*  Do the Wick contractions of the quark propagators to form the 
 *  two point functions
 *
 *  This routine assumes that the HQET propagatotors have had their 
 *  source gamma matrix representation flipped.
 *
 *   Function arguments
 * 
 *    On input
 *        q_zonked      :: site structure pointer to spin_wilson_vector 
 *                         [[ zonked quark ]]
 *        q_sequential  :: site structure pouinter to spin_wilson_vector 
 *                         [[ sequential inverted quark ]]
 *        zonked_pt :: which zonked quark 
 *        spect_pt  :: which spectator quark
 *
 *    On output
 *        corr :: the two point function correlators 
 *
 *  
 *   This code calls the general purpose function  "meson_cont_mom",
 *   which does the low level contraction. This routine includes the 
 *   hermitian conjugate and multiplies the zonked quark by gamma_5
 *   on both dides.
 */


void contract_light_twopt(complex *corr, field_offset q_zonked, 
			  field_offset q_sequential,
			  int zonked_pt, int spect_pt)
{
  double t_start ;
  int base_pt, q_stride, op_stride;

  t_start = dclock() ;

  /* Compute partial offset for storage of result in corr[] */

  base_pt   = TWOPT_FORM_WHERE(0,zonked_pt,spect_pt,0,0 )  ; 
  q_stride  = TWOPT_FORM_WHERE(0,0,        0,       1,0 )  ; 
  op_stride = TWOPT_FORM_WHERE(0,0,        0,       0,1 )  ; 

  meson_cont_mom(corr , q_zonked, q_sequential, 
		 base_pt, q_stride, op_stride,
		 two_pt, MAX_TWOPT);

  IF_VERBOSE_ON(1)
    printf("Time to Wick contract light 2pt correlators = %g sec\n",
	   dclock() - t_start) ;


} 



