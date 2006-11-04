/******* load_qop_asqtad_coeffs.c ****************/
/* MIMD version 7 */

/* Loads coefficient structure for QOP routines */

/* CD 10/07/06 Split from fermion_force_asqtad_qop.c */

#include "generic_ks_includes.h"
#include <qop.h>

void load_qop_asqtad_coeffs(QOP_asqtad_coeffs_t *c, Real weight,
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

