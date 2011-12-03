/******* load_qop_asqtad_coeffs_P.c ****************/
/* MIMD version 7 */

/* Loads coefficient structure for QOP routines */

/* Note: This is an include file for load_qop_asqtad_coeffs_F.c and
   load_qop_asqtad_coeffs_D.c */

/* CD 10/07/06 Split from fermion_force_asqtad_qop.c */

#if ( QOP_Precision == 1 )
#define LOAD_QOP_ASQTAD_COEFFS  load_qop_F_asqtad_coeffs
#else
#define LOAD_QOP_ASQTAD_COEFFS  load_qop_D_asqtad_coeffs
#endif

#include "generic_ks_includes.h"
#include "../include/generic_ks_qop.h"
#include <qop.h>

void LOAD_QOP_ASQTAD_COEFFS(QOP_asqtad_coeffs_t *c, Real weight,
			    Real *act_path_coeff)
{
  Real ferm_epsilon;

  ferm_epsilon = 2.0*weight;
  
  /* Path coefficients times fermion epsilon */

  c->one_link     = act_path_coeff[0]*ferm_epsilon ; 
  c->naik         = act_path_coeff[1]*ferm_epsilon ;
  c->three_staple = act_path_coeff[2]*ferm_epsilon ;
  c->five_staple  = act_path_coeff[3]*ferm_epsilon ;
  c->seven_staple = act_path_coeff[4]*ferm_epsilon ;
  c->lepage       = act_path_coeff[5]*ferm_epsilon ;
}

