/******* ks_action_coeffs_hisq_qop.c ****************/
/* MIMD version 7 */

/* Loads HISQ coefficient structure for QOP routines */

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
#include "../include/umethod.h"
#include <qop.h>

/*--------------------------------------------------------------------*/

static QOP_hisq_unitarize_group_t 
convert_milc_to_qop_ugroup(int ug_milc){
  QOP_hisq_unitarize_group_t ug_qop = QOP_UNITARIZE_U3;
  char myname[] = "convert_milc_to_qop_ugroup";

  switch(ug_milc){
  case UNITARIZE_U3:
    ug_qop = QOP_UNITARIZE_U3;
    break;
  case UNITARIZE_SU3:
    ug_qop = QOP_UNITARIZE_SU3;
    break;
  default:
    if(mynode()==0)
      printf("%s: Unknown unitarization group. Set UNITARIZE_U3 or UNITARIZE_SU3\n",
	     myname);
    terminate(1);
  }
  return ug_qop;
}

static QOP_hisq_unitarize_method_t
convert_milc_to_qop_umethod(int um_milc){
  QOP_hisq_unitarize_method_t um_qop = UNITARIZE_ANALYTIC;
  char myname[] = "convert_milc_to_qop_umethod";

  /* (Add cases when they are supported in QOP) */

  switch(um_milc){
    //  case UNITARIZE_NONE:
    //    um_qop = QOP_UNITARIZE_NONE;
    //    break;
    //  case UNITARIZE_ROOT:
    //    um_qop = QOP_UNITARIZE_ROOT;
    //    break;
  case UNITARIZE_RATIONAL:
  case UNITARIZE_ANALYTIC:
    um_qop = QOP_UNITARIZE_RATIONAL;
    break;
  default:
    if(mynode()==0)printf("%s: Unknown unitarization method %d.\n", myname, um_milc);
    terminate(1);
  }
  return um_qop;
}

QOP_hisq_coeffs_t *
create_hisq_coeffs_qop(ks_action_paths_hisq *ap){

  int i;
  QOP_hisq_coeffs_t *ac;
  int n_naiks = get_n_naiks(ap);
  double *eps_naik = get_eps_naik(ap);
  int ugroup = get_ugroup(ap);
  int umethod = get_umethod(ap);
  
  char myname[] = "create_hisq_coeffs_qop";

  ac = (QOP_hisq_coeffs_t *)malloc(sizeof(QOP_hisq_coeffs_t));
  if(ac == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  if(n_naiks > QOP_MAX_NAIK){
    printf("%s: n_naiks = %d exceeds QOP_MAX_NAIK = %d\n",
	   myname, n_naiks, QOP_MAX_NAIK);
    terminate(1);
  }

  ac->n_naiks = n_naiks;
  for(i = 0; i < n_naiks; i++)
    ac->eps_naik[i] = eps_naik[i];

  /* WARNING! We are setting the weight to 1 here.  Then the links are
     created with the correct conventions for the inverter.  But for
     the fermion force calculation, any flavor weights, such as
     2*nflavor/4 must be included in the fermion epsilon factors when
     calling the fermion force term.
  */

  ac->fat7_one_link       = ap->p1.act_path_coeff.one_link ;
  ac->fat7_three_staple   = ap->p1.act_path_coeff.three_staple ;
  ac->fat7_five_staple    = ap->p1.act_path_coeff.five_staple ;
  ac->fat7_seven_staple   = ap->p1.act_path_coeff.seven_staple ;

  ac->asqtad_one_link     = ap->p2.act_path_coeff.one_link ;
  ac->asqtad_naik         = ap->p2.act_path_coeff.naik ;
  ac->asqtad_three_staple = ap->p2.act_path_coeff.three_staple ;
  ac->asqtad_five_staple  = ap->p2.act_path_coeff.five_staple ;
  ac->asqtad_seven_staple = ap->p2.act_path_coeff.seven_staple ;
  ac->asqtad_lepage       = ap->p2.act_path_coeff.lepage ;

  ac->difference_one_link = ap->p3.act_path_coeff.one_link ;
  ac->difference_naik     = ap->p3.act_path_coeff.naik ;

  ac->ugroup  = convert_milc_to_qop_ugroup(ugroup);
  ac->umethod = convert_milc_to_qop_umethod(umethod);

  return ac;
}

void 
destroy_hisq_coeffs_qop(QOP_hisq_coeffs_t *ac){
  if(ac == NULL)return;
  free(ac);
}


/*--------------------------------------------------------------------*/

int get_n_naiks_qop(QOP_hisq_coeffs_t *ac){
  return ac->n_naiks;
}

double *get_eps_naik_qop(QOP_hisq_coeffs_t *ac){
  return ac->eps_naik;
}

int get_ugroup_qop(QOP_hisq_coeffs_t *ac){
  return ac->ugroup;
}

#define MAX_STRING 2048

char *
get_ac_string_hisq_qop(QOP_hisq_coeffs_t *ac){
  static char str[MAX_STRING] = "";
  
  snprintf(str, MAX_STRING, "action.hisq.fat7.one_link %e\naction.hisq.fat7.three_staple %e\naction.hisq.fat7.five_staple %e\naction.hisq.fat7.seven_staple %e\naction.hisq.asqtad.one_link %e\naction.hisq.asqtad.three_staple %e\naction.hisq.asqtad.five_staple %e\naction.hisq.asqtad.seven_staple %e\naction.hisq.asqtad.lepage %e\naction.hisq.asqtad.naik %e\naction.hisq.difference.one_link %e\naction.hisq.difference.naik %e\n",
	   ac->fat7_one_link, ac->fat7_three_staple, ac->fat7_five_staple, ac->fat7_seven_staple, 
	   ac->asqtad_one_link, ac->asqtad_three_staple, ac->asqtad_five_staple, ac->asqtad_seven_staple, 
	   ac->asqtad_lepage, ac->asqtad_naik, ac->difference_one_link, ac->difference_naik);

  str[MAX_STRING-1] = '\0';
  return str;
}

