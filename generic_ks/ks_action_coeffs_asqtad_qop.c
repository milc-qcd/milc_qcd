/******* ks_action_coeffs_asqtad_qop.c ****************/
/* MIMD version 7 */

/* Loads asqtad coefficient structure for QOP routines */

/* Note: This is an include file for load_qop_asqtad_coeffs_F.c and
   load_qop_asqtad_coeffs_D.c */

/* CD 10/07/06 Split from fermion_force_asqtad_qop.c */
/* CD 4/11 reworked from load_qop_asqtad_coeffs_qop_P.c
   and load_qop_hisq_coeffs_qop_P.c */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../include/comdefs.h"
#include "../include/ks_action_paths.h"
#include "../include/ks_action_coeffs_qop.h"
#include <qop.h>

/*--------------------------------------------------------------------*/

QOP_asqtad_coeffs_t *
create_asqtad_coeffs_qop(ks_action_paths *ap){

  QOP_asqtad_coeffs_t *ac;
  asqtad_coeffs_t *act_path_coeff;
  char myname[] = "create_asqtad_coeffs_qop";

  if(ap == NULL)return NULL;

  act_path_coeff = &ap->p.act_path_coeff;

  ac = (QOP_asqtad_coeffs_t *)malloc(sizeof(QOP_asqtad_coeffs_t));
  if(ac == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  /* WARNING! We are setting the weight to 1 here.  Then the links are
     created with the correct conventions for the inverter.  But for
     the fermion force calculation, any flavor weights, such as
     2*nflavor/4 must be included in the fermion epsilon factors when
     calling the fermion force term.
  */

  ac->one_link     = act_path_coeff->one_link     ;
  ac->naik         = act_path_coeff->naik         ;
  ac->three_staple = act_path_coeff->three_staple ;
  ac->five_staple  = act_path_coeff->five_staple  ;
  ac->seven_staple = act_path_coeff->seven_staple ;
  ac->lepage       = act_path_coeff->lepage       ;

  return ac;
}

void 
destroy_asqtad_coeffs_qop(QOP_asqtad_coeffs_t *ac){
  if(ac == NULL)return;
  free(ac);
}

#define MAX_STRING 2048

char *
get_ac_string_asqtad_qop(QOP_asqtad_coeffs_t *ac){
  static char str[MAX_STRING] = "";

  snprintf(str, MAX_STRING, "action.asqtad.one_link %e\naction.asqtad.three_staple %e\naction.asqtad.five_staple %e\naction.asqtad.seven_staple %e\naction.asqtad.lepage %e\naction.asqtad.naik %e\n",
	   ac->one_link, ac->three_staple, ac->five_staple, ac->seven_staple, ac->lepage, ac->naik);

  str[MAX_STRING-1] = '\0';
  return str;
}
