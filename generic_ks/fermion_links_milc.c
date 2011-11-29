/****** fermion_links_milc.c  ******************/
/* MIMD version 7 */

/* Mar 2011 CD */

#include <stdlib.h>
#include <stdio.h>
#include "../include/comdefs.h"
#include "../include/fermion_links.h"
#include "../include/fermion_links_milc.h"
#include "../include/ks_action_paths.h"
#include "../include/info.h"

/*--------------------------------------------------*/
/* Create/destroy the fm_ap_links_t structure      */
/*--------------------------------------------------*/

static fm_ap_links_t *
create_fm_ap_links_t(info_t *info, ks_action_paths *ap, su3_matrix *links, 
		     ferm_links_options_t *options){

  fm_ap_links_t *al;
  char myname[] = "create_fm_ap_links_t";

  al = (fm_ap_links_t *)malloc(sizeof(fm_ap_links_t));
  if(al == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }
  
  al->ap = ap;
  al->fm = create_imp_ferm_links();
  load_imp_ferm_links(info, al->fm, ap, links, options->want_back);

  return al;
}

static void
destroy_fm_ap_links_t(fm_ap_links_t *al){
  if(al == NULL)return;

  destroy_imp_ferm_links(al->fm);
  destroy_path_table(al->ap);

  free(al);
}

static void
invalidate_fm_ap_links_t(fm_ap_links_t *al){
  if(al == NULL)return;

  destroy_imp_ferm_links(al->fm);
  al->fm = NULL;
}

static void
restore_fm_ap_links_t(info_t *info, fm_ap_links_t *al, su3_matrix *links, 
		      int want_back){
  if(al == NULL)return;

  /* If the fm member is not NULL, the links are assumed to be
     valid */

  if(al->fm != NULL)return;

  al->fm = create_imp_ferm_links();
  load_imp_ferm_links(info, al->fm, al->ap, links, want_back);
}

static imp_ferm_links_t **
get_fm_ap_links_t_fm(fm_ap_links_t *al){
  return &al->fm;
}

static ks_action_paths*
get_fm_ap_links_ap(fm_ap_links_t *al){
  return al->ap;
}

/* We keep only one precision for MILC types */
static int
valid_fm_ap_links_t(fm_ap_links_t *al, int precision){
  return al->fm != NULL;
}

/*-------------------------------------------------------------------*/
/* Create/destroy the milc_fm_links_t structure */
/*-------------------------------------------------------------------*/

static milc_fm_links_t *
create_milc_fm_links_t(info_t *info, ks_action_paths *ap, ks_action_paths *ap_du0, 
		       su3_matrix *links, ferm_links_options_t *options ){

  milc_fm_links_t *al;
  char myname[] = "create_milc_fm_links_t";
  double final_flop = 0.0;

  al = (milc_fm_links_t *)malloc(sizeof(milc_fm_links_t));
  if(al == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  al->fm_ap = create_fm_ap_links_t(info, ap, links, options);
  final_flop += info->final_flop;

  if(options->want_du0){
    al->fm_ap_du0 = create_fm_ap_links_t(info, ap_du0, links, options);
    final_flop += info->final_flop;
  } else {
    al->fm_ap_du0 = NULL;
  }

  info->final_flop = final_flop;
  return al;
}

static void 
destroy_milc_fm_links_t(milc_fm_links_t *al){
  if(al == NULL) return;

  destroy_fm_ap_links_t(al->fm_ap);
  destroy_fm_ap_links_t(al->fm_ap_du0);

  free(al);
}

static void
invalidate_milc_fm_links_t(milc_fm_links_t *al){
  if(al == NULL)return;
  invalidate_fm_ap_links_t(al->fm_ap);
  invalidate_fm_ap_links_t(al->fm_ap_du0);
}

static void
restore_milc_fm_links_t(info_t *info, milc_fm_links_t *al, 
			    su3_matrix *links, int want_back){
  double final_flop = 0.0;
  info->final_flop = 0.0;
  if(al == NULL)return;

  restore_fm_ap_links_t(info, al->fm_ap, links, want_back);
  final_flop += info->final_flop;

  restore_fm_ap_links_t(info, al->fm_ap_du0, links, want_back);
  final_flop += info->final_flop;

  info->final_flop = final_flop;
}

static imp_ferm_links_t **
get_milc_fm_ap_links_fm(milc_fm_links_t *al){
  return get_fm_ap_links_t_fm(al->fm_ap);
}

static imp_ferm_links_t **
get_milc_fm_ap_du0_links_fm(milc_fm_links_t *al){
  return get_fm_ap_links_t_fm(al->fm_ap_du0);
}

static ks_action_paths *
get_milc_fm_ap_links_ap(milc_fm_links_t *al){
  return get_fm_ap_links_ap(al->fm_ap);
}

static int
valid_milc_fm_ap_links(milc_fm_links_t *al, int precision){
  /* We assume that if the fm_ap links are valid,
     so are the fm_ap_du0 links */
  return valid_fm_ap_links_t(al->fm_ap, precision);
}

/*********************************************************************/
/* The public API                                                    */
/*********************************************************************/

/*-----------------------------------------------------------------*/
/* Create the fermion_links structure for improved fermion actions */
/*-----------------------------------------------------------------*/

#ifdef FLTIME
static const char *milc_prec[2] = {"F", "D"};
#endif

fermion_links_t *
create_fermion_links(int precision, int phases_in, su3_matrix *links){
  
  fermion_links_t *fl;
  ks_action_paths *ap, *ap_du0;
  info_t info = INFO_ZERO;
  char myname[] = "create_fermion_links";

  /* Precision for MILC is ignored: use the prevailing precision */

  if(precision != PRECISION)
    if(mynode() == 0)printf("%s: Warning. Precision request replaced by %d\n",
			    myname, PRECISION);

  if( phases_in != 1){
    if(mynode() == 0)printf("BOTCH: %s needs phases in\n",myname); 
    terminate(1);
  }
  
  fl = create_fermion_links_t();

  /* Create the path tables */

  /* (We copy the pointers into the fm_ap_links_t objects
     and the responsibility for freeing space is handed over to
     "destroy_fm_ap_links_t") */

  ap = create_path_table();
  if(fl->options.want_du0)
    ap_du0 = create_path_table();
  else
    ap_du0 = NULL;

  make_path_table(ap, ap_du0);

  /* Complete the structure */

  fl->flg = create_milc_fm_links_t(&info, ap, ap_du0, links, &fl->options);


#ifdef FLTIME
  if(mynode()==0)printf("FLTIME: time = %e (asqtad %s) mflops = %e\n",
	       info.final_sec,milc_prec[PRECISION-1],
	       info.final_flop/(1e6*info.final_sec) );
#endif
  return fl;
}

/*----------------------------------------*/
/* Destroy the fermion links structure */
/*----------------------------------------*/

void 
destroy_fermion_links(fermion_links_t *fl){

  if(fl == NULL)return;

  destroy_milc_fm_links_t(fl->flg);
  destroy_fermion_links_t(fl);
}

/*----------------------------------------*/
/* Invalidate the fermion links structure */
/*----------------------------------------*/

void 
invalidate_fermion_links(fermion_links_t *fl){
  if(fl == NULL)return;
  invalidate_milc_fm_links_t(fl->flg);
}

/*----------------------------------------*/
/* Restore links that may be invalid      */
/*----------------------------------------*/

void 
restore_fermion_links(fermion_links_t *fl, int precision, int phases_in, su3_matrix *links){

  char myname[] = "restore_fermion_links_fm";
  info_t info = INFO_ZERO;

  if(fl == NULL){
    if(mynode() == 0)printf("%s: Called with a null pointer\n", myname);
    terminate(1);
  }

  if(precision != PRECISION)
    if(mynode() == 0)printf("%s: Warning. Precision request replaced by %d\n",
			    myname, PRECISION);

  if( phases_in != 1){
    if(mynode() == 0)printf("BOTCH: %s needs phases in\n",myname); terminate(1);
  }
  
  restore_milc_fm_links_t(&info, fl->flg, links, fl->options.want_back);
#ifdef FLTIME
  if(mynode()==0)printf("FLTIME: time = %e (asqtad %s) mflops = %e\n",
	       info.final_sec,milc_prec[PRECISION-1],
	       info.final_flop/(1e6*info.final_sec) );
#endif
}

/*----------------------------------------*/
/* Accessors                              */
/*----------------------------------------*/

imp_ferm_links_t **
get_fm_links(fermion_links_t *fl){
  return get_milc_fm_ap_links_fm(fl->flg);
}

imp_ferm_links_t **
get_fm_du0_links(fermion_links_t *fl){
  return get_milc_fm_ap_du0_links_fm(fl->flg);
}

ks_action_paths *
get_action_paths(fermion_links_t *fl){
  return get_milc_fm_ap_links_ap(fl->flg);
}

int 
valid_fermion_links(fermion_links_t *fl, int precision){
  if(fl == NULL)return 0;
  return valid_milc_fm_ap_links(fl->flg, precision);
}

char *
get_action_parameter_string(fermion_links_t *fl){
  ks_action_paths *ap = get_action_paths(fl);
  char *str = get_ap_string(ap);

  return str;
}
