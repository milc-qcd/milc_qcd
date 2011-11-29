/****** fermion_links_asqtad_qop.c  *************************/
/* MIMD version 7 */

/* Apr 2011 CD */

#include <stdlib.h>
#include <stdio.h>
#include <qop.h>
#include "../include/comdefs.h"
#include "../include/fermion_links.h"
#include "../include/fermion_links_qop.h"
#include "../include/ks_action_coeffs_qop.h"
#include "../include/generic_qop.h"

/*--------------------------------------------------*/
/* Create/destroy the qop_fm_ac_links_t structure      */
/*--------------------------------------------------*/

static qop_fm_ac_links_t *
create_qop_fm_ac_links_t(QOP_asqtad_coeffs_t *ac, int precision,
			 su3_matrix *links, ferm_links_options_t *options){

  qop_fm_ac_links_t *al;
  char myname[] = "create_qop_fm_ac_links_t";

  al = (qop_fm_ac_links_t *)malloc(sizeof(qop_fm_ac_links_t));
  if(al == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }
  
  al->ac = ac;
  al->fm = create_fn_links_qop();
  load_fn_links_qop(al->fm, al->ac, precision, links, options->want_back);

  return al;
}

static void
destroy_qop_fm_ac_links_t(qop_fm_ac_links_t *al){
  if(al == NULL)return;

  destroy_fn_links_qop(al->fm);
  destroy_asqtad_coeffs_qop(al->ac);

  free(al);
}

static void
invalidate_qop_fm_ac_links_t(qop_fm_ac_links_t *al){
  if(al == NULL)return;

  destroy_fn_links_qop(al->fm);
  al->fm = NULL;
}

static void
restore_qop_fm_ac_links_t(qop_fm_ac_links_t *al, int precision, 
			  su3_matrix *links, int want_back){
  if(al == NULL)return;

  /* If we already have links of the requested precision, those links
     are assumed to be valid, and we don't recompute them.  If we have
     links of the wrong precision, we destroy them and create new
     ones. */

  if(al->fm != NULL){
    if((precision == 1 && get_F_asqtad_links(al->fm) != NULL) ||
       (precision == 2 && get_D_asqtad_links(al->fm) != NULL) )return;
  }

  destroy_fn_links_qop(al->fm);
  al->fm = create_fn_links_qop();
  load_fn_links_qop(al->fm, al->ac, precision, links, want_back);
}

static imp_ferm_links_t **
get_qop_fm_ac_links_t_fm(qop_fm_ac_links_t *al){
  return &al->fm;
}

static QOP_asqtad_coeffs_t *
get_qop_fm_ac_links_ac(qop_fm_ac_links_t *al){
  return al->ac;
}

static int
valid_qop_fm_ac_links_t(qop_fm_ac_links_t *al, int precision){
  if(al->fm != NULL){
    if((precision == 1 && get_F_asqtad_links(al->fm) != NULL) ||
       (precision == 2 && get_D_asqtad_links(al->fm) != NULL) )
      return 1;
  }

  return 0;
}

/*-------------------------------------------------------------------*/
/* Create/destroy the qop_asqtad_links_t structure */
/*-------------------------------------------------------------------*/

static qop_asqtad_links_t *
create_qop_asqtad_links_t(QOP_asqtad_coeffs_t *ac, QOP_asqtad_coeffs_t *ac_du0, 
			  int precision, su3_matrix *links, 
			  ferm_links_options_t *options ){

  qop_asqtad_links_t *al;
  char myname[] = "create_qop_asqtad_links_t";

  al = (qop_asqtad_links_t *)malloc(sizeof(qop_asqtad_links_t));
  if(al == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  al->fm_ac = create_qop_fm_ac_links_t(ac, precision, links, options);

  if(options->want_du0)
    al->fm_ac_du0 = create_qop_fm_ac_links_t(ac_du0, precision, links, options);
  else
    al->fm_ac_du0 = NULL;

  return al;
}

static void 
destroy_qop_asqtad_links_t(qop_asqtad_links_t *al){
  if(al == NULL) return;

  destroy_qop_fm_ac_links_t(al->fm_ac);
  destroy_qop_fm_ac_links_t(al->fm_ac_du0);

  free(al);
}

static void
invalidate_qop_asqtad_links_t(qop_asqtad_links_t *al){
  if(al == NULL)return;
  invalidate_qop_fm_ac_links_t(al->fm_ac);
  invalidate_qop_fm_ac_links_t(al->fm_ac_du0);
}

static void
restore_qop_asqtad_links_t(qop_asqtad_links_t *al, int precision,
			    su3_matrix *links, int want_back){
  if(al == NULL)return;

  restore_qop_fm_ac_links_t(al->fm_ac, precision, links, want_back);
  restore_qop_fm_ac_links_t(al->fm_ac_du0, precision, links, want_back);
}

static imp_ferm_links_t **
get_qop_asqtad_ac_links_fm(qop_asqtad_links_t *al){
  if(al == NULL)return NULL;
  return get_qop_fm_ac_links_t_fm(al->fm_ac);
}

static imp_ferm_links_t **
get_qop_asqtad_ac_du0_links_fm(qop_asqtad_links_t *al){
  if(al == NULL)return NULL;
  return get_qop_fm_ac_links_t_fm(al->fm_ac_du0);
}

static QOP_asqtad_coeffs_t *
get_qop_asqtad_ac_links_ac(qop_asqtad_links_t *al){
  if(al == NULL)return NULL;
  return get_qop_fm_ac_links_ac(al->fm_ac);
}

static int
valid_qop_asqtad_ac_links(qop_asqtad_links_t *al, int precision){
  /* We assume that if the qop_fm_ac links are valid,
     so are the qop_fm_ac_du0 links */
  return valid_qop_fm_ac_links_t(al->fm_ac, precision);
}

/*********************************************************************/
/* The public API                                                    */
/*********************************************************************/

/*-----------------------------------------------------------------*/
/* Create the fermion_links structure for improved fermion actions */
/*-----------------------------------------------------------------*/

fermion_links_t *
create_fermion_links(int precision, int phases_in, su3_matrix *links){
  
  fermion_links_t *fl;
  QOP_asqtad_coeffs_t *ac, *ac_du0;
  ks_action_paths *ap, *ap_du0;
  char myname[] = "create_fermion_links";
  

  if( phases_in != 1){
    if(mynode() == 0)printf("BOTCH: %s needs phases in\n",myname); 
    terminate(1);
  }

  /* Make sure QOP is initialized */

  if(initialize_qop() != QOP_SUCCESS){
    printf("%s: Failed to initialize qop\n",myname);
    terminate(1);
  }

  /* Create the parent structure */

  fl = create_fermion_links_t();

  /* (We copy the pointers into the fn_ap_links_t objects
     and the responsibility for freeing space is handed over to
     "destroy_hisq_links_t") */

  /* (We copy the pointers into the fn_ap_links_t objects
     and the responsibility for freeing space is handed over to
     "destroy_hisq_links_t") */

  /* Create the path tables */

  ap = create_path_table();
  if(fl->options.want_du0)
    ap_du0 = create_path_table();
  else
    ap_du0 = NULL;

  /* TO DO: don't need the paths here -- just the coeffs.
     So distinguish setting coeffs and creating paths. */
  make_path_table(ap, ap_du0);

  /* Use them to create the action coefficients */

  ac = create_asqtad_coeffs_qop(ap);
  ac_du0 = create_asqtad_coeffs_qop(ap_du0);

  /* Destroy the MILC path tables */

  destroy_path_table(ap);
  destroy_path_table(ap_du0);

  /* Complete the structure */

  fl->flg = create_qop_asqtad_links_t(ac, ac_du0, precision, links, &fl->options);

  return fl;
}

/*----------------------------------------*/
/* Destroy the fermion links structure */
/*----------------------------------------*/

void 
destroy_fermion_links(fermion_links_t *fl){

  if(fl == NULL)return;

  destroy_qop_asqtad_links_t(fl->flg);
  destroy_fermion_links_t(fl);
}

/*----------------------------------------*/
/* Invalidate the fermion links structure */
/*----------------------------------------*/

void 
invalidate_fermion_links(fermion_links_t *fl){
  if(fl == NULL)return;
  invalidate_qop_asqtad_links_t(fl->flg);
}

/*----------------------------------------*/
/* Restore links that may be invalid      */
/*----------------------------------------*/

void 
restore_fermion_links(fermion_links_t *fl, int precision,
		      int phases_in, su3_matrix *links){

  char myname[] = "restore_fermion_links_fm";

  if(fl == NULL){
    if(mynode() == 0)printf("%s: Called with a null pointer\n", myname);
    terminate(1);
  }

  if( phases_in != 1){
    if(mynode() == 0)printf("BOTCH: %s needs phases in\n",myname); terminate(1);
  }
  
  restore_qop_asqtad_links_t(fl->flg, precision, links, fl->options.want_back);
}

/*----------------------------------------*/
/* Accessors                              */
/*----------------------------------------*/

imp_ferm_links_t **
get_fm_links(fermion_links_t *fl){
  return get_qop_asqtad_ac_links_fm(fl->flg);
}

imp_ferm_links_t **
get_fm_du0_links(fermion_links_t *fl){
  return get_qop_asqtad_ac_du0_links_fm(fl->flg);
}

QOP_asqtad_coeffs_t *
get_action_coeffs(fermion_links_t *fl){
  return get_qop_asqtad_ac_links_ac(fl->flg);
}

int 
valid_fermion_links(fermion_links_t *fl, int precision){
  if(fl == NULL)return 0;
  return valid_qop_asqtad_ac_links(fl->flg, precision);
}

char *
get_action_parameter_string(fermion_links_t *fl){
  QOP_asqtad_coeffs_t *ac = get_action_coeffs(fl);
  char *str = get_ac_string_asqtad_qop(ac);

  return str;
}

/********************************************************************/
