/******************** contract_two_pt_prop.c ************************/
/* MIMD version 6 */
/*   TWO POINT FUNCTIONS FOR PROPAGATING FORM FACTORS
 *
 *  A separate function for the heavy-light and light-light
 *  correlators.
 */  


/*  
 *  Do the Wick contractions for the relevant two point functions
 *
 *  This routine will need to be updated, when the light quark
 *  smearing is implemented properly
 */


/* Modifications

   C. DeTar 4/21/98 (lean) Provision for storing values only for time slices on node

 */

#include "prop_form_includes.h"
#include "opertypes.h"
#include "corrlist.h"

#ifdef DEBUGDEF
#include "debug_form.h"
#endif


/*              LIGHT---LIGHT routine 
 *  Do the Wick contractions of the quark propagators to form the 
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
 *        spect_pt  :: which spactator quark
 *
 *    On output
 *        corr :: the two point function correlators 
 *
 *  
 *   This code calls the general purpose function  "meson_cont_mom_lean2",
 *   which does the low level contraction. This routine includes the 
 *   hermitian conjugate and multiplies the zonked quark by gamma_5
 *   on both sides.
 *
 *   corr = trace(\Gamma_in \gamma_5 q_zonked^{\dagger} \gamma_5 \Gamma_out q_spectator) 
 *
 */


void contract_LL2(complex *corr, field_offset q_zonked, field_offset q_spectator,
int zonked_pt, int spect_pt)
{
  int base_pt, q_stride, op_stride;
  double t_start ;

  t_start = dclock() ;


  /************************************************************/

  /* Compute partial offset for storage of result in corr[] */

  base_pt   = LL_TWOPT_FORM_WHERE(0,zonked_pt,spect_pt,0,0 )  ; 
  q_stride  = LL_TWOPT_FORM_WHERE(0,0,        0,       1,0 )  ; 
  op_stride = LL_TWOPT_FORM_WHERE(0,0,        0,       0,1 )  ; 

  meson_cont_mom_lean2(corr , q_zonked, q_spectator,
		 base_pt, q_stride, op_stride,
		 w_meson_store_t,w_meson_my_t,w_meson_nstore,
		 no_k_values,k_momstore,
		 MAX_TWOPT, two_pt,
		 F_OFFSET(QTMP),DIMQTMP);

  IF_VERBOSE_ON(1)
    printf("Time to Wick contract light-light 2pt correlators = %g sec\n",
	   dclock() - t_start) ;


} 



/*
 *  Contraction for the heavy-light two point functions.
 *
 *  The spectator quark propagator is contracted with the 
 *  heavy zonked quark propagator.
 *
 *   corr = trace(\Gamma_in \gamma_5 q_zonked^{\dagger} \gamma_5 \Gamma_out q_spectator) 
 *
 */


void contract_HL2(complex *corr, field_offset q_zonked, field_offset q_spectator,
int zonked_pt, int spect_pt)
{
  int base_pt, q_stride, op_stride;
  double t_start ;

  t_start = dclock() ;


  /************************************************************/

  /* Compute partial offset for storage of result in corr[] */

  base_pt   =    HL_TWOPT_FORM_WHERE(0,zonked_pt,spect_pt,0,0 )  ; 
  q_stride  =    HL_TWOPT_FORM_WHERE(0,0,        0,       1,0 )  ; 
  op_stride =    HL_TWOPT_FORM_WHERE(0,0,        0,       0,1 )  ; 

  meson_cont_mom_lean2(corr , q_zonked, q_spectator,
		 base_pt, q_stride, op_stride,
		 w_meson_store_t,w_meson_my_t,w_meson_nstore,
		 no_k_values,k_momstore,
		 MAX_TWOPT, two_pt,
		 F_OFFSET(QTMP),DIMQTMP);

  IF_VERBOSE_ON(1)
    printf("Time to Wick contract heavy-light 2pt correlators = %g sec\n",
	   dclock() - t_start) ;


} 









/*
 *  Contraction for the heavy-light two point functions with operator
 *  insertions.
 *
 *  The spectator quark propagator is contracted with the 
 *  heavy zonked quark propagator.
 *
 *   corr = trace(\Gamma_in \gamma_5 q_zonked^{\dagger} \gamma_5 \Gamma_out q_spectator) 
 *
 */


void contract_HL2_with_rotations(complex *corr, field_offset q_zonked, field_offset q_spectator,
field_offset q_rot, int zonked_pt, int spect_pt)
{
  int base_pt, q_stride, op_stride;
  double t_start ;

  t_start = dclock() ;


  /************************************************************/

  /* Compute partial offset for storage of result in corr[] */

  base_pt   =    HL_TWOPT_FORM_WHERE(0,zonked_pt,spect_pt,0,0 )  ; 
  q_stride  =    HL_TWOPT_FORM_WHERE(0,0,        0,       1,0 )  ; 
  op_stride =    HL_TWOPT_FORM_WHERE(0,0,        0,       0,1 )  ; 

  meson_cont_mom_lean2(corr , q_zonked, q_spectator,
		 base_pt, q_stride, op_stride,
		 w_meson_store_t,w_meson_my_t,w_meson_nstore,
		 no_k_values,k_momstore,
		 MAX_TWOPT, two_pt,
		 F_OFFSET(QTMP),DIMQTMP);


  /*** 3D derivative on the spectator ***/
  base_pt += op_stride*MAX_TWOPT;
  clover_rotate_fermilab(q_spectator, q_rot);

  meson_cont_mom_lean2(corr , q_zonked, q_rot,
		 base_pt, q_stride, op_stride,
		 w_meson_store_t,w_meson_my_t,w_meson_nstore,
		 no_k_values,k_momstore,
		 MAX_TWOPT, two_pt,
		 F_OFFSET(QTMP),DIMQTMP);

  /*** 3D derivative on the zonked ***/
  base_pt += op_stride*MAX_TWOPT;
  clover_rotate_fermilab(q_zonked, q_rot);

  meson_cont_mom_lean2(corr , q_rot, q_spectator,
		 base_pt, q_stride, op_stride,
		 w_meson_store_t,w_meson_my_t,w_meson_nstore,
		 no_k_values,k_momstore,
		 MAX_TWOPT, two_pt,
		 F_OFFSET(QTMP),DIMQTMP);


  /*** 4D derivative on the spectator ***/
  /** base_pt += op_stride*MAX_TWOPT;
  clover_rotate(q_spectator, q_rot);

  meson_cont_mom_lean2(corr , q_zonked, q_rot,
		 base_pt, q_stride, op_stride,
		 w_meson_store_t,w_meson_my_t,w_meson_nstore,
		 no_k_values,k_momstore,
		 MAX_TWOPT, two_pt,
		 F_OFFSET(QTMP),DIMQTMP); **/


  /*** 4D derivative on the zonked ***/
  /** base_pt += op_stride*MAX_TWOPT;
  clover_rotate(q_zonked, q_rot);

  meson_cont_mom_lean2(corr , q_rot, q_spectator,
		 base_pt, q_stride, op_stride,
		 w_meson_store_t,w_meson_my_t,w_meson_nstore,
		 no_k_values,k_momstore,
		 MAX_TWOPT, two_pt,
		 F_OFFSET(QTMP),DIMQTMP); **/

  /************************************************************/

  IF_VERBOSE_ON(1)
    printf("Time to Wick contract heavy-light 2pt correlators + operators = %g sec\n",
	   dclock() - t_start) ;

} 





