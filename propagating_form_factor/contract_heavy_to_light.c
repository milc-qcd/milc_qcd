/********************** contract_heavy_to_light.c **********************/
/* MIMD version 6 */
/*
 *  Do the Wick contractions for the heavy-to-light
 *  three point function.
 */

/* Modifications

   C. DeTar 4/30/97 Now working from a table of source-sink gamma matrix pairs
   C. DeTar 4/21/98 (lean) Provision for storing values only for time slices on node

 */

#include "prop_form_includes.h"
#include "opertypes.h"
#include "corrlist.h"

#ifdef DEBUGDEF
#include "debug_form.h"
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
 *        q_sequential  :: site structure pouinter to spin_wilson_vector 
 *                         [[ sequential inverted quark ]]
 *        zonked_pt :: which zonked quark 
 *        spact_pt  :: which spactator quark
 *        p_pt      :: momentum of the B meson
 *        seq_pt    :: the sequential source counter
 *
 *    On output
 *        corr :: the three point function correlators 
 *
 *
 *   This code calls the general purpose function  "meson_cont_mom_lean",
 *   which does the low level contraction. This routine includes the 
 *   hermitian conjugate and multiplies the zonked quark by gamma_5
 *   on both dides.
 *
 */


void contract_HL3(complex *corr, field_offset q_zonked, field_offset q_sequential,
field_offset q_rot, int p_pt, int seq_pt, int zonked_pt, int spect_pt)
{
  double t_start ;
  int base_pt, q_stride, op_stride ;

  /**********----------**********----------**********/


  t_start = dclock() ;

  /* Compute partial offset for storage of result in corr[] */
 
  base_pt    = LIGHT_FORM_WHERE(0,zonked_pt,seq_pt,spect_pt,0,p_pt, 0 ) ; 
  q_stride   = LIGHT_FORM_WHERE(0,0,0,        0,       1,0,      0 ) ;
  op_stride  = LIGHT_FORM_WHERE(0,0,0,        0,       0,0,      1 ) ;

  /* First, contract zonked and sequential */

  meson_cont_mom_lean2(corr , q_zonked, q_sequential, 
		 base_pt, q_stride, op_stride, 
		 w_meson_store_t,w_meson_my_t,w_meson_nstore,
		 no_q_values,q_momstore,
		 MAX_THREEPT, three_pt,
		 F_OFFSET(QTMP),DIMQTMP);



  /*** 3D derivative on the sequential ***/
  base_pt += op_stride*MAX_THREEPT;
  clover_rotate_fermilab(q_sequential, q_rot);

  meson_cont_mom_lean2(corr, q_zonked, q_rot, 
		 base_pt, q_stride, op_stride, 
		 w_meson_store_t,w_meson_my_t,w_meson_nstore,
		 no_q_values,q_momstore,
		 MAX_THREEPT, three_pt,
		 F_OFFSET(QTMP),DIMQTMP);

  /*** 3D derivative on the zonked ***/
  base_pt += op_stride*MAX_THREEPT;
  clover_rotate_fermilab(q_zonked, q_rot);

  meson_cont_mom_lean2(corr, q_rot, q_sequential, 
		 base_pt, q_stride, op_stride, 
		 w_meson_store_t,w_meson_my_t,w_meson_nstore,
		 no_q_values,q_momstore,
		 MAX_THREEPT, three_pt,
		 F_OFFSET(QTMP),DIMQTMP);



  /*** 4D derivative on the sequential ***/
  /**  base_pt += op_stride*MAX_THREEPT;
  clover_rotate(q_sequential, q_rot);

  meson_cont_mom_lean2(corr, q_zonked, q_rot, 
		 base_pt, q_stride, op_stride, 
		 w_meson_store_t,w_meson_my_t,w_meson_nstore,
		 no_q_values,q_momstore,
		 MAX_THREEPT, three_pt,
		 F_OFFSET(QTMP),DIMQTMP); **/

  /*** 4D derivative on the zonked ***/
  /** base_pt += op_stride*MAX_THREEPT;
  clover_rotate(q_zonked, q_rot);

  meson_cont_mom_lean2(corr, q_rot, q_sequential, 
		 base_pt, q_stride, op_stride, 
		 w_meson_store_t,w_meson_my_t,w_meson_nstore,
		 no_q_values,q_momstore,
		 MAX_THREEPT, three_pt,
		 F_OFFSET(QTMP),DIMQTMP); **/




  /**** end of the section ****/

  IF_VERBOSE_ON(1)
    printf("Time to Wick contract all the heavy-light 3pt correlators = %g sec\n",
	   dclock() - t_start) ;

}
