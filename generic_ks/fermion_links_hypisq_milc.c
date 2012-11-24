/****** fermion_links_hypisq_milc.c  ******************/

/* PLACEHOLDER COPIED FROM fermion_links_hisq_milc.c NEEDS REVISION
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
/* Create/destroy the hypisq_links_t structure      */
/*--------------------------------------------------*/
 
static hypisq_links_t *
create_hypisq_links_t(info_t *info, ks_action_paths_hypisq *ap, su3_matrix *links, 
		      ferm_links_options_t *options){

  /* Allocate the structure and create the auxiliary links and fn
     links, based on the given path tables and gauge links.
  */

  hypisq_links_t *hl;
  int i;
  char myname[] = "create_hypisq_links_t";

  hl = (hypisq_links_t *)malloc(sizeof(hypisq_links_t));
  if(hl == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  for(i = 0; i < MAX_NAIK; i++)
    hl->fn[i] = NULL;
  
  hl->ap = ap;

  create_hypisq_links_milc(info, hl->fn, &hl->fn_deps, &hl->aux, ap, links, 
			   options->want_deps, options->want_back);

  /* Free the space if so desired */
  if(!options->want_aux){
    destroy_hypisq_auxiliary_t(hl->aux);
    hl->aux = NULL;
  }

  return hl;
}

static void
destroy_hypisq_links_t(hypisq_links_t *hl){
  
  if(hl == NULL) return;

  /* Destroy all members and free the structure */

  destroy_hypisq_links_milc(hl->ap, hl->aux, hl->fn, hl->fn_deps);
  destroy_path_table_hypisq(hl->ap);

  free(hl);
}


static void
invalidate_hypisq_links_t(hypisq_links_t *hl){
  if(hl == NULL)return;

  /* Destroy the auxiliary and fn links and reset the pointers */
  /* Keep the path tables */

  destroy_hypisq_links_milc(hl->ap, hl->aux, hl->fn, hl->fn_deps);
  hl->aux = NULL;
  hl->fn_deps = NULL;
}

static void
restore_hypisq_links_t(info_t *info, hypisq_links_t *hl, su3_matrix *links, 
		       ferm_links_options_t *options){

  if(hl == NULL)return;

  /* If the first fn member is not NULL, the links are assumed to be
     valid and we do nothing */

  if(hl->fn[0] != NULL)return;

  /* Allocate and create the HYPISQ auxiliary links and the fn links */

  create_hypisq_links_milc(info, hl->fn, &hl->fn_deps, &hl->aux, 
			   hl->ap, links, options->want_deps,
			   options->want_back);

  /* Free the space if so desired */
  if(!options->want_aux){
    destroy_hypisq_auxiliary_t(hl->aux);
    hl->aux = NULL;
  }
}

static fn_links_t **
get_hypisq_links_t_fn(hypisq_links_t *hl){
  if(hl == NULL)return NULL;
  return hl->fn;
}

static fn_links_t *
get_hypisq_links_t_fn_deps(hypisq_links_t *hl){
  if(hl == NULL)return NULL;
  return hl->fn_deps;
}

static ks_action_paths_hypisq*
get_hypisq_links_t_ap(hypisq_links_t *hl){
  if(hl == NULL)return NULL;
  return hl->ap;
}

/* We keep only one precision for MILC types */
static int
valid_hypisq_links_t(hypisq_links_t *hl, int precision){
  return hl->fn[0] != NULL;
}

static hypisq_auxiliary_t *
get_hypisq_links_t_aux(hypisq_links_t *hl){
  if(hl == NULL)return NULL;
  return hl->aux;
}


/*-------------------------------------------------------------------*/
/* Create/destroy the milc_hypisq_links_t structure */
/*-------------------------------------------------------------------*/

static milc_hypisq_links_t *
create_milc_hypisq_links_t(info_t *info, ks_action_paths_hypisq *ap, 
			   su3_matrix *links, 
			   ferm_links_options_t *options){

  milc_hypisq_links_t *hl;
  char myname[] = "create_milc_hypisq_links_t";

  hl = (milc_hypisq_links_t *)malloc(sizeof(milc_hypisq_links_t));
  if(hl == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  hl->hypisq = create_hypisq_links_t(info, ap, links, options);

  return hl;
}

static void 
destroy_milc_hypisq_links_t(milc_hypisq_links_t *hl){
  if(hl == NULL)return;

  destroy_hypisq_links_t(hl->hypisq);

  free(hl);
}

static void
invalidate_milc_hypisq_links_t(milc_hypisq_links_t *hl){
  if(hl == NULL)return;
  invalidate_hypisq_links_t(hl->hypisq);
}

static void
restore_milc_hypisq_links_t(info_t *info, milc_hypisq_links_t *hl, 
			    su3_matrix *links, 
			    ferm_links_options_t *options){
  if(hl == NULL)return;

  restore_hypisq_links_t(info, hl->hypisq, links, options);
}

static fn_links_t **
get_milc_hypisq_links_fn(milc_hypisq_links_t *hl){
  if(hl == NULL)return NULL;
  return get_hypisq_links_t_fn(hl->hypisq);
}

static fn_links_t *
get_milc_hypisq_links_fn_deps(milc_hypisq_links_t *hl){
  if(hl == NULL)return NULL;
  return get_hypisq_links_t_fn_deps(hl->hypisq);
}

static ks_action_paths_hypisq *
get_milc_hypisq_links_ap(milc_hypisq_links_t *hl){
  if(hl == NULL)return NULL;
  return get_hypisq_links_t_ap(hl->hypisq);
}

static int
valid_milc_hypisq_links(milc_hypisq_links_t *hl, int precision){
  /* We assume that if the hypisq links are valid,
     so are the hypisq_du0 links */
  if(hl == NULL)return 0;
  return valid_hypisq_links_t(hl->hypisq, precision);
}

static hypisq_auxiliary_t *
get_milc_hypisq_links_aux(milc_hypisq_links_t *hl){
  if(hl == NULL)return NULL;
  return get_hypisq_links_t_aux(hl->hypisq);
}

/*********************************************************************/
/* The public API                                                    */
/*********************************************************************/

/*----------------------------------------------------*/
/* Create the fermion_links structure for hypisq      */
/*----------------------------------------------------*/

#ifdef FLTIME
static const char *milc_prec[2] = {"F", "D"};
#endif

fermion_links_t *
create_fermion_links_hypisq(int precision, int n_naiks, 
			    double eps_naik[], int phases_in, su3_matrix *links){
  
  fermion_links_t *fl;
  ks_action_paths_hypisq *ap;
  info_t info = INFO_ZERO;
  char myname[] = "create_fermion_links_hypisq";

  /* Precision for MILC is ignored: use the prevailing precision */

  if(precision != PRECISION)
    if(mynode() == 0)printf("%s: Warning. Precision request replaced by %d\n", myname,
		 PRECISION);

  if( phases_in != 1){
    if(mynode() == 0)printf("BOTCH: %s needs phases in\n",myname); terminate(1);
  }
  
  fl = create_fermion_links_t();

  /* Create the path tables */

  /* (We copy the pointers into the fn_ap_links_t objects
     and the responsibility for freeing space is handed over to
     "destroy_hypisq_links_t") */

  ap = create_path_table_hypisq();

  make_path_table_hypisq(ap, n_naiks, eps_naik);

  /* Complete the structure */

  fl->flg = create_milc_hypisq_links_t(&info, ap, links, &fl->options);

#ifdef FLTIME
  if(mynode()==0)printf("FLTIME: time = %e (HYPISQ %s) mflops = %e\n",
			info.final_sec,milc_prec[PRECISION-1],
			info.final_flop/(1e6*info.final_sec) );
#endif
  return fl;
}

/*----------------------------------------*/
/* Destroy the fermion links structure */
/*----------------------------------------*/

void 
destroy_fermion_links_hypisq(fermion_links_t *fl){

  if(fl == NULL) return;

  destroy_milc_hypisq_links_t(fl->flg);
  destroy_fermion_links_t(fl);
}


/*----------------------------------------*/
/* Invalidate the fermion links structure */
/*----------------------------------------*/

void 
invalidate_fermion_links(fermion_links_t *fl){
  if(fl == NULL)return;
  invalidate_milc_hypisq_links_t(fl->flg);
}

/*----------------------------------------*/
/* Restore links that may be invalid      */
/*----------------------------------------*/

/* If the links are valid, do nothing.  Otherwise,
   restore them from the gauge links provided. 
   The KS and periodic/antiperiodic boundary
   phases must be in the links.

*/

void 
restore_fermion_links_hypisq(fermion_links_t *fl, int precision, 
			     int phases_in, su3_matrix *links){

  char myname[] = "restore_fermion_links";
  info_t info = INFO_ZERO;

  if(fl == NULL){
    if(mynode() == 0)printf("%s: Called with a null pointer\n", myname);
    terminate(1);
  }

  if(precision != PRECISION)
    if(mynode() == 0)printf("%s: Warning. Precision request replaced by %d\n",
			    myname, PRECISION);

  if( phases_in != 1){
    if(mynode() == 0)printf("BOTCH: %s needs phases in\n",myname); 
    terminate(1);
  }
  
  restore_milc_hypisq_links_t(&info, fl->flg, links, &fl->options);

#ifdef FLTIME
  if(mynode()==0)printf("FLTIME: time = %e (HYPISQ %s) mflops = %e\n",
			info.final_sec,milc_prec[PRECISION-1],
			info.final_flop/(1e6*info.final_sec) );
#endif
}

/*----------------------------------------*/
/* Accessors                              */
/*----------------------------------------*/

fn_links_t **
get_fm_links(fermion_links_t *fl){
  if(fl == NULL)return NULL;
  return get_milc_hypisq_links_fn(fl->flg);
}

fn_links_t *
get_fn_deps_links(fermion_links_t *fl){
  if(fl == NULL)return NULL;
  return get_milc_hypisq_links_fn_deps(fl->flg);
}

ks_action_paths_hypisq *
get_action_paths_hypisq(fermion_links_t *fl){
  if(fl == NULL)return NULL;
  return get_milc_hypisq_links_ap(fl->flg);
}

int
get_n_naiks_hypisq(fermion_links_t *fl){
  ks_action_paths_hypisq* ap;
  if(fl == NULL)return -1;
  ap = get_action_paths_hypisq(fl);
  return get_n_naiks(ap);
}

double *
get_eps_naik_hypisq(fermion_links_t *fl){
  ks_action_paths_hypisq* ap;
  if(fl == NULL)return NULL;
  ap = get_action_paths_hypisq(fl);
  return get_eps_naik(ap);
}

int 
valid_fermion_links(fermion_links_t *fl, int precision){
  if(fl == NULL)return 0;
  return valid_milc_hypisq_links(fl->flg, precision);
}

hypisq_auxiliary_t *
get_hypisq_auxiliary(fermion_links_t *fl){
  if(fl == NULL)return NULL;
  return get_milc_hypisq_links_aux(fl->flg);
}

char *
get_action_parameter_string(fermion_links_t *fl){
  ks_action_paths_hypisq *ap = get_action_paths_hypisq(fl);
  char *str = get_ap_string_hypisq(ap);

  return str;
}
