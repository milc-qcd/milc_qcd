/****** fermion_links_hisq_milc.c  ******************/
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
/* Create/destroy the hisq_links_t structure      */
/*--------------------------------------------------*/
 
static hisq_links_t *
create_hisq_links_t(info_t *info, ks_action_paths_hisq *ap, su3_matrix *links, 
		    ferm_links_options_t *options){

  /* Allocate the structure and create the auxiliary links and fn
     links, based on the given path tables and gauge links.
  */

  hisq_links_t *hl;
  int i;
  char myname[] = "create_hisq_links_t";

  hl = (hisq_links_t *)malloc(sizeof(hisq_links_t));
  if(hl == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  hl->fn0 = NULL;
  hl->fn_deps = NULL;
  hl->ap = ap;

  create_hisq_links_milc(info, &hl->fn0, &hl->fn_deps, &hl->aux, ap, links, 
			 options->want_deps, options->want_back);

  /* Free the space if so desired */
  if(!options->want_aux){
    destroy_hisq_auxiliary_t(hl->aux);
    hl->aux = NULL;
  }

  return hl;
}

static void
destroy_hisq_links_t(hisq_links_t *hisq){
  
  if(hisq == NULL) return;

  /* Destroy all members and free the structure */

  destroy_hisq_links_milc(hisq->ap, hisq->aux, hisq->fn0, hisq->fn_deps);
  destroy_path_table_hisq(hisq->ap);
  hisq->aux = NULL;
  hisq->fn0 = NULL;
  hisq->fn_deps = NULL;
}


static void
invalidate_hisq_links_t(hisq_links_t *hisq){
  if(hisq == NULL)return;

  /* Destroy the auxiliary and fn links and reset the pointers */
  /* Keep the path tables */

  destroy_hisq_links_milc(hisq->ap, hisq->aux, hisq->fn0, hisq->fn_deps);
  hisq->aux = NULL;
  hisq->fn0 = NULL;
  hisq->fn_deps = NULL;
}

static void
restore_hisq_links_t(info_t *info, hisq_links_t *hisq, su3_matrix *links, 
		     ferm_links_options_t *options){

  if(hisq == NULL)return;

  /* If the fn0 and aux members ares not NULL, we assume the links are
     valid and we do nothing */

  int valid = 0;
  if(options->want_aux)
    valid = hisq->aux != NULL && hisq->fn0 != NULL;
  else
    valid = hisq->fn0 != NULL;
  if(valid)return;

  /* Allocate and create the HISQ auxiliary links and the fn links */

  create_hisq_links_milc(info, &hisq->fn0, &hisq->fn_deps, &hisq->aux, 
			 hisq->ap, links, options->want_deps,
			 options->want_back);

  /* Free the space if so desired */
  if(!options->want_aux){
    destroy_hisq_auxiliary_t(hisq->aux);
    hisq->aux = NULL;
  }
}

static ks_action_paths_hisq*
get_hisq_links_t_ap(hisq_links_t *hl){
  if(hl == NULL)return NULL;
  return hl->ap;
}

/* We keep only one precision for MILC types */
static int
valid_hisq_links_t(hisq_links_t *hl, int precision){
  return hl->fn0 != NULL;
}

/* We keep only one precision for MILC types */
static int
valid_hisq_aux_t(hisq_links_t *hisq, int precision){
  return hisq->aux != NULL;
}

static hisq_auxiliary_t *
get_hisq_links_t_aux(hisq_links_t *hisq){
  if(hisq == NULL)return NULL;
  return hisq->aux;
}


/*-------------------------------------------------------------------*/
/* Create/destroy the milc_hisq_links_t structure */
/*-------------------------------------------------------------------*/

static milc_hisq_links_t *
create_milc_hisq_links_t(info_t *info, ks_action_paths_hisq *ap, 
			 su3_matrix *links, 
			 ferm_links_options_t *options){

  milc_hisq_links_t *hl;
  char myname[] = "create_milc_hisq_links_t";

  hl = (milc_hisq_links_t *)malloc(sizeof(milc_hisq_links_t));
  if(hl == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  hl->hisq = create_hisq_links_t(info, ap, links, options);

  return hl;
}

static void 
destroy_milc_hisq_links_t(milc_hisq_links_t *flg){
  if(flg == NULL)return;

  destroy_hisq_links_t(flg->hisq);
}

static void
invalidate_milc_hisq_links_t(milc_hisq_links_t *flg){
  if(flg == NULL)return;
  invalidate_hisq_links_t(flg->hisq);
}

static void
restore_milc_hisq_links_t(info_t *info, milc_hisq_links_t *flg, 
			  su3_matrix *links, 
			  ferm_links_options_t *options){
  if(flg == NULL)return;

  restore_hisq_links_t(info, flg->hisq, links, options);
}

static ks_action_paths_hisq *
get_milc_hisq_links_ap(milc_hisq_links_t *flg){
  if(flg == NULL)return NULL;
  return get_hisq_links_t_ap(flg->hisq);
}

static int
valid_milc_hisq_links(milc_hisq_links_t *flg, int precision){
  /* We assume that if the hisq links are valid,
     so are the hisq_du0 links */
  if(flg == NULL)return 0;
  return valid_hisq_links_t(flg->hisq, precision);
}

static int
valid_milc_hisq_aux(milc_hisq_links_t *flg, int precision){
  if(flg->hisq == NULL)return 0;
  else  return valid_hisq_aux_t(flg->hisq, precision);
}

static fn_links_t *
get_hisq_links_t_fn_deps(hisq_links_t *hisq){
  if(hisq == NULL)return NULL;
  return hisq->fn_deps;
}

static hisq_auxiliary_t *
get_milc_hisq_links_aux(milc_hisq_links_t *hl){
  if(hl == NULL)return NULL;
  return get_hisq_links_t_aux(hl->hisq);
}

static fn_links_t *
get_milc_hisq_links_fn_deps(milc_hisq_links_t *hl){
  if(hl == NULL)return NULL;
  return get_hisq_links_t_fn_deps(hl->hisq);
}

int
get_n_naiks_hisq(fermion_links_t *fl){
  ks_action_paths_hisq* ap;
  if(fl == NULL)return -1;
  ap = get_action_paths_hisq(fl);
  return get_n_naiks(ap);
}

static fn_links_t *
get_hisq_links_t_fn(hisq_links_t *hl, int i_naik, ferm_links_options_t *options){
  if(hl == NULL)return NULL;
  int n_naiks = hl->ap->n_naiks;
  if(i_naik >= n_naiks){
    node0_printf("ERROR: Requested Naik epsilon index %d >= n_naik = %d\n", i_naik, n_naiks);
    terminate(1);
  }

  fn_links_t *fn;

  fn_links_t *fn_deps = get_hisq_links_t_fn_deps(hl);
  double *eps_naik = hl->ap->eps_naik;

  if(i_naik == 0){

    fn = hl->fn0;
    fn->preserve = 1;
  } else {
    fn = create_fn_links();
    if(options->want_back){
      fn->fatback = create_fatlinks();
      fn->lngback = create_lnglinks();
    }
    scalar_mult_fn(fn_deps, eps_naik[i_naik], fn);
    add_fn(fn, hl->fn0, fn);
    fn->preserve = 0;
    fn->eps_naik = eps_naik[i_naik];
  }
  return fn;
}

static fn_links_t *
get_milc_hisq_links_fn(milc_hisq_links_t *hl, int i_naik, ferm_links_options_t *options){
  if(hl == NULL)return NULL;
  return get_hisq_links_t_fn(hl->hisq, i_naik, options);
}

/*********************************************************************/
/* The public API                                                    */
/*********************************************************************/

/*----------------------------------------------------*/
/* Create the fermion_links structure for hisq      */
/*----------------------------------------------------*/

#ifdef FLTIME
static const char *milc_prec[2] = {"F", "D"};
#endif

fermion_links_t *
create_fermion_links_hisq(int precision, int n_naiks, 
			  double eps_naik[], int phases_in, su3_matrix *links){
  
  fermion_links_t *fl;
  ks_action_paths_hisq *ap;
  info_t info = INFO_ZERO;
  char myname[] = "create_fermion_links_hisq";

  /* Precision for MILC is ignored: use the prevailing precision */

  if(precision != MILC_PRECISION)
    if(mynode() == 0)printf("%s: Warning. Precision request replaced by %d\n", myname,
		 MILC_PRECISION);

  if( phases_in != 1){
    if(mynode() == 0) printf("BOTCH: %s needs phases in\n",myname);
    terminate(1);
  }
  
  fl = create_fermion_links_t();

  /* Create the path tables */

  /* (We copy the pointers into the fn_ap_links_t objects
     and the responsibility for freeing space is handed over to
     "destroy_hisq_links_t") */

  ap = create_path_table_hisq();

  make_path_table_hisq(ap, n_naiks, eps_naik);

  /* Complete the structure */

  fl->flg = create_milc_hisq_links_t(&info, ap, links, &fl->options);

#ifdef FLTIME
#ifdef USE_FL_GPU

#if defined(HAVE_QUDA)
  if(mynode()==0)printf("FLTIME: time = %e (HISQ QUDA %s) mflops = %e\n",
	       info.final_sec,milc_prec[MILC_PRECISION-1],
	       info.final_flop/(1e6*info.final_sec) );
#elif defined(HAVE_GRID)
  if(mynode()==0)printf("FLTIME: time = %e (HISQ GRID %s) mflops = %e\n",
	       info.final_sec,milc_prec[MILC_PRECISION-1],
	       info.final_flop/(1e6*info.final_sec) );
#endif

#else
  if(mynode()==0)printf("FLTIME: time = %e (HISQ MILC %s) mflops = %e\n",
	       info.final_sec,milc_prec[MILC_PRECISION-1],
	       info.final_flop/(1e6*info.final_sec) );

#endif
#endif
  return fl;
}

/*----------------------------------------*/
/* Destroy the fermion links structure */
/*----------------------------------------*/

void 
destroy_fermion_links_hisq(fermion_links_t *fl){

  if(fl == NULL) return;

  destroy_milc_hisq_links_t(fl->flg);
  destroy_fermion_links_t(fl);
}


/*----------------------------------------*/
/* Invalidate the fermion links structure */
/*----------------------------------------*/

void invalidate_fermion_links(fermion_links_t *fl){
  if(fl == NULL)return;
  invalidate_milc_hisq_links_t(fl->flg);
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
restore_fermion_links_hisq(fermion_links_t *fl, int precision, 
			   int phases_in, su3_matrix *links){

  char myname[] = "restore_fermion_links";
  info_t info = INFO_ZERO;

  if(fl == NULL){
    if(mynode() == 0)printf("%s: Called with a null pointer\n", myname);
    terminate(1);
  }

  if(precision != MILC_PRECISION)
    if(mynode() == 0)printf("%s: Warning. Precision request replaced by %d\n",
			    myname, MILC_PRECISION);

  if( phases_in != 1){
    if(mynode() == 0)printf("BOTCH: %s needs phases in\n",myname); 
    terminate(1);
  }
  
  restore_milc_hisq_links_t(&info, fl->flg, links, &fl->options);

#ifdef FLTIME
#ifdef USE_FL_GPU
  if(mynode()==0)printf("FLTIME: time = %e (HISQ QUDA %s) mflops = %e\n",
	       info.final_sec,milc_prec[MILC_PRECISION-1],
	       info.final_flop/(1e6*info.final_sec) );
#else
  if(mynode()==0)printf("FLTIME: time = %e (HISQ MILC %s) mflops = %e\n",
	       info.final_sec,milc_prec[MILC_PRECISION-1],
	       info.final_flop/(1e6*info.final_sec) );
#endif
#endif
}

/*----------------------------------------*/
/* Accessors                              */
/*----------------------------------------*/

fn_links_t *
get_fm_links(fermion_links_t *fl, int i_naik){
  if(fl == NULL)return NULL;
  ferm_links_options_t options = fl->options;
  return get_milc_hisq_links_fn(fl->flg, i_naik, &options);
}

fn_links_t *
get_fn_deps_links(fermion_links_t *fl){
  if(fl == NULL)return NULL;
  return get_milc_hisq_links_fn_deps(fl->flg);
}

ks_action_paths_hisq *
get_action_paths_hisq(fermion_links_t *fl){
  if(fl == NULL)return NULL;
  return get_milc_hisq_links_ap(fl->flg);
}

double *
get_eps_naik_hisq(fermion_links_t *fl){
  ks_action_paths_hisq* ap;
  if(fl == NULL)return NULL;
  ap = get_action_paths_hisq(fl);
  return get_eps_naik(ap);
}

int 
valid_fermion_links(fermion_links_t *fl, int precision){
  if(fl == NULL)return 0;
  int valid = 0;
  if(fl->options.want_aux == 1)
    valid = valid_milc_hisq_aux(fl->flg, precision) && valid_milc_hisq_links(fl->flg, precision);
  else
    valid =  valid_milc_hisq_links(fl->flg, precision);
  return valid;
}

hisq_auxiliary_t *
get_hisq_auxiliary(fermion_links_t *fl){
  if(fl == NULL)return NULL;
  return get_milc_hisq_links_aux(fl->flg);
}

char *
get_action_parameter_string(fermion_links_t *fl){
  ks_action_paths_hisq *ap = get_action_paths_hisq(fl);
  char *str = get_ap_string_hisq(ap);

  return str;
}
