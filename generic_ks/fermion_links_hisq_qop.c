/****************** fermion_links_hisq_qop.c ***********************/
/* MIMD version 7 */

#include <stdlib.h>
#include <stdio.h>
#include "../include/comdefs.h"
#include "../include/fermion_links.h"
#include "../include/fermion_links_qop.h"
#include "../include/fn_links.h"
#include "../include/ks_action_paths.h"
#include "../include/ks_action_coeffs_qop.h"
#include "../include/generic_qop.h"

/*---------------------------------------------------*/
/* Create/destroy the qop_hisq_ac_links_t structure  */
/*---------------------------------------------------*/
 
static qop_hisq_ac_links_t *
create_qop_hisq_ac_links_t(QOP_hisq_coeffs_t *ac, int precision,
			   su3_matrix *links, ferm_links_options_t *options){

  int i;
  qop_hisq_ac_links_t *hisq;
  char myname[] = "create_qop_hisq_ac_links_t";

  hisq = (qop_hisq_ac_links_t *)malloc(sizeof(qop_hisq_ac_links_t));
  if(hisq == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }
  
  hisq->ac = ac;
  hisq->hl = create_hisq_links_qop(ac, precision, links, options->want_deps,
				   options->want_aux);

  /* The fn_links here are needed for MILC wrappers.
     They are structures with only lists of copied pointers. */

  for(i = 0; i < ac->n_naiks; i++){
    hisq->fn[i] = create_fn_links_qop();
    set_asqtad_links_from_hisq(hisq->fn[i], hisq->hl, i);
  }

  if(options->want_deps){
    hisq->fn_deps = create_fn_links_qop();
    set_asqtad_deps_links_from_hisq(hisq->fn_deps, hisq->hl);
  }
  else
    hisq->fn_deps = NULL;

  return hisq;
}

/* Unset copied pointers and free allocated fn link data */
static void
unset_qop_hisq_fn_links(qop_hisq_ac_links_t *hisq){
  int i;
  for(i = 0; i < hisq->ac->n_naiks; i++){
    /* (We must first clear copied pointers to prevent double-freeing
       allocated memory) */
    unset_asqtad_links_from_hisq(hisq->fn[i]);
  }

  unset_asqtad_deps_links_from_hisq(hisq->fn_deps);
}

static void
destroy_qop_hisq_ac_links_t(qop_hisq_ac_links_t *hisq){

  int i;
  
  if(hisq == NULL) return;

  /* Destroy all members and free the structure */

  for(i = 0; i < hisq->ac->n_naiks; i++){
    /* (We must first clear copied pointers to prevent double-freeing
       allocated memory) */
    unset_asqtad_links_from_hisq(hisq->fn[i]);
    destroy_fn_links_qop(hisq->fn[i]);
  }

  unset_asqtad_deps_links_from_hisq(hisq->fn_deps);
  destroy_fn_links_qop(hisq->fn_deps);
  hisq->fn_deps = NULL;

  destroy_hisq_links_qop(hisq->hl);

  free(hisq);
}


static void
invalidate_qop_hisq_ac_links_t(qop_hisq_ac_links_t *hisq){
  int i;

  if(hisq == NULL)return;

  /* Destroy the hisq links, but keep the action coefficients */

  destroy_hisq_links_qop(hisq->hl);
  hisq->hl = NULL;

  /* Clear copied pointers but don't free the parent structures */

  for(i = 0; i < hisq->ac->n_naiks; i++)
    unset_asqtad_links_from_hisq(hisq->fn[i]);

  unset_asqtad_deps_links_from_hisq(hisq->fn_deps);
  
}

static void
restore_qop_hisq_ac_links_t(qop_hisq_ac_links_t *hisq, int precision,
			    su3_matrix *links, ferm_links_options_t *options){
  int i;
  if(hisq == NULL)return;

  /* If we already have links of the requested precision, those links
     are assumed to be valid, and we don't recompute them.  If we have
     links of the wrong precision, we destroy them and create new
     ones. */

  if(hisq->hl != NULL){
    if((precision == 1 && get_F_hisq_links_qop(hisq->hl) != NULL) ||
       (precision == 2 && get_D_hisq_links_qop(hisq->hl) != NULL) )return;
  }

  /* Allocate and create the HISQ auxiliary links and the fn links */

  unset_qop_hisq_fn_links(hisq);
  destroy_hisq_links_qop(hisq->hl);
  hisq->hl = create_hisq_links_qop(hisq->ac, precision, links, options->want_deps,
				   options->want_aux);

  for(i = 0; i < hisq->ac->n_naiks; i++)
    set_asqtad_links_from_hisq(hisq->fn[i], hisq->hl, i);

  set_asqtad_deps_links_from_hisq(hisq->fn_deps, hisq->hl);

}

static fn_links_qop_t **
get_qop_hisq_ac_links_t_fn(qop_hisq_ac_links_t *hisq){

  if(hisq == NULL)return NULL;
  if(hisq->hl == NULL)return NULL;

  return hisq->fn;
}

static fn_links_qop_t *
get_qop_hisq_ac_links_t_fn_deps(qop_hisq_ac_links_t *hisq){
  if(hisq == NULL)return NULL;
  if(hisq->hl == NULL)return NULL;
  if(hisq->fn_deps == NULL)return NULL;

  return hisq->fn_deps;
}

static QOP_F3_FermionLinksHisq *
get_F_qop_hisq_ac_links_t_hl(qop_hisq_ac_links_t *hisq){
  if(hisq == NULL)return NULL;
  return get_F_hisq_links_qop(hisq->hl);
}

static QOP_D3_FermionLinksHisq *
get_D_qop_hisq_ac_links_t_hl(qop_hisq_ac_links_t *hisq){
  if(hisq == NULL)return NULL;
  return get_D_hisq_links_qop(hisq->hl);
}

static QOP_hisq_coeffs_t*
get_qop_hisq_ac_links_t_ac(qop_hisq_ac_links_t *hisq){
  if(hisq == NULL)return NULL;
  return hisq->ac;
}

/* For now we don't distinguish precisions */
static int
valid_qop_hisq_ac_links_t(qop_hisq_ac_links_t *hisq, int precision){
  if(hisq->hl != NULL){
    if((precision == 1 && get_F_hisq_links_qop(hisq->hl) != NULL) ||
       (precision == 2 && get_D_hisq_links_qop(hisq->hl) != NULL) )
      return 1;
  }
  return 0;
}

/*-------------------------------------------------------------------*/
/* Create/destroy the qop_hisq_links_t structure                     */
/*-------------------------------------------------------------------*/

static qop_hisq_links_t *
create_qop_hisq_links_t(QOP_hisq_coeffs_t *ac, 
			int precision, su3_matrix *links, 
			ferm_links_options_t *options){

  qop_hisq_links_t *hisq;
  char myname[] = "create_qop_hisq_links_t";

  hisq = (qop_hisq_links_t *)malloc(sizeof(qop_hisq_links_t));
  if(hisq == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  hisq->hisq_ac = create_qop_hisq_ac_links_t(ac, precision, links, options);

  return hisq;
}

static void 
destroy_qop_hisq_links_t(qop_hisq_links_t *hisq){
  if(hisq == NULL)return;

  destroy_qop_hisq_ac_links_t(hisq->hisq_ac);

  free(hisq);
}

static void
invalidate_qop_hisq_links_t(qop_hisq_links_t *hisq){
  if(hisq == NULL)return;

  invalidate_qop_hisq_ac_links_t(hisq->hisq_ac);
}

static void
restore_qop_hisq_links_t(qop_hisq_links_t *hisq, 
			 int precision, su3_matrix *links, 
			 ferm_links_options_t *options){
  if(hisq == NULL)return;

  restore_qop_hisq_ac_links_t(hisq->hisq_ac, precision, links, options);
}

static fn_links_qop_t **
get_qop_hisq_links_fn(qop_hisq_links_t *hisq){
  if(hisq == NULL)return NULL;
  return get_qop_hisq_ac_links_t_fn(hisq->hisq_ac);
}

static fn_links_qop_t *
get_qop_hisq_links_fn_deps(qop_hisq_links_t *hisq){
  if(hisq == NULL)return NULL;
  return get_qop_hisq_ac_links_t_fn_deps(hisq->hisq_ac);
}

QOP_F3_FermionLinksHisq *
get_F_qop_hisq_links_hl(qop_hisq_links_t *hisq){
  if(hisq == NULL)return NULL;
  return get_F_qop_hisq_ac_links_t_hl(hisq->hisq_ac);
}

QOP_D3_FermionLinksHisq *
get_D_qop_hisq_links_hl(qop_hisq_links_t *hisq){
  if(hisq == NULL)return NULL;
  return get_D_qop_hisq_ac_links_t_hl(hisq->hisq_ac);
}

static QOP_hisq_coeffs_t *
get_qop_hisq_links_ac(qop_hisq_links_t *hisq){
  if(hisq == NULL)return NULL;
  return get_qop_hisq_ac_links_t_ac(hisq->hisq_ac);
}

static int
valid_qop_hisq_links(qop_hisq_links_t *hisq, int precision){
  if(hisq == NULL)return 0;
  return valid_qop_hisq_ac_links_t(hisq->hisq_ac, precision);
}

static int
set_qop_hisq_link_options(void){
  int status;

  /* Note: the want_deps and want_aux options are set in hisq_links_qop.c */

  /* Set values */
  QOP_opt_t qop_hl_opt[5] = {

#ifdef HISQ_REUNIT_ALLOW_SVD
    {.tag = "reunit_allow_svd",.value=1},
#else
    {.tag = "reunit_allow_svd",.value=0},
#endif

#ifdef HISQ_REUNIT_SVD_ONLY
    {.tag = "reunit_svd_only",.value=1},
#else
    {.tag = "reunit_svd_only",.value=0},
#endif

#ifdef HISQ_REUNIT_SVD_REL_ERROR
    {.tag = "reunit_svd_rel_error",.value=HISQ_REUNIT_SVD_REL_ERROR},
#else
    {.tag = "reunit_svd_rel_error",.value=1e-8},
#endif

#ifdef HISQ_REUNIT_SVD_ABS_ERROR
    {.tag = "reunit_svd_abs_error",.value=HISQ_REUNIT_SVD_ABS_ERROR},
#else
    {.tag = "reunit_svd_abs_error",.value=1e-8},
#endif

#ifdef HISQ_SVD_VALUES_INFO
    {.tag = "svd_values_info",.value=1}
#else
    {.tag = "svd_values_info",.value=0}
#endif

  };

  /* Set options */

  status = QOP_hisq_links_set_opts(qop_hl_opt, 5);

  if(status!= QOP_SUCCESS)
    if(mynode()==0)printf("set_qop_hisq_link_options: ERROR setting QOP options\n");

  return status;
}

/*********************************************************************/
/* The public API                                                    */
/*********************************************************************/

/*----------------------------------------------------*/
/* Create the fermion_links structure for hisq      */
/*----------------------------------------------------*/

fermion_links_t *
create_fermion_links_hisq(int precision, int n_naiks, 
			  double eps_naik[], int phases_in, 
			  su3_matrix *links){
  
  fermion_links_t *fl;
  QOP_hisq_coeffs_t *ac;
  ks_action_paths_hisq *ap;
  char myname[] = "create_fermion_links_hisq";

  if( phases_in != 1){
    if(mynode() == 0)printf("BOTCH: %s needs phases in\n",myname); 
    terminate(1);
  }
  
  /* Make sure QOP is initialized */

  if(initialize_qop() != QOP_SUCCESS){
    printf("%s: Failed to initialize qop\n",myname);
    terminate(1);
  }

  /* Set QOP HISQ link options */

  set_qop_hisq_link_options();

  /* Create the parent structure */

  fl = create_fermion_links_t();

  /* Create the path table */

  /* (We copy the pointers into the fn_ap_links_t objects
     and the responsibility for freeing space is handed over to
     "destroy_hisq_links_t") */

  ap = create_path_table_hisq();

  /* Get the action coefficients, but don't construct the paths */
  /* (QOP will construct them if needed) */

  load_act_path_coeff_hisq(ap, n_naiks, eps_naik);

  /* Use them to create the action coefficients */

  ac = create_hisq_coeffs_qop(ap);

  /* Destroy the MILC path tables */

  destroy_path_table_hisq(ap);

  /* Complete the structure */

  fl->flg = create_qop_hisq_links_t(ac, precision, links, &fl->options);

  return fl;
}

/*----------------------------------------*/
/* Destroy the fermion links structure */
/*----------------------------------------*/

void 
destroy_fermion_links_hisq(fermion_links_t *fl){

  if(fl == NULL) return;

  destroy_qop_hisq_links_t(fl->flg);
  destroy_fermion_links_t(fl);
}


/*----------------------------------------*/
/* Invalidate the fermion links structure */
/*----------------------------------------*/

void 
invalidate_fermion_links(fermion_links_t *fl){
  if(fl == NULL)return;
  invalidate_qop_hisq_links_t(fl->flg);
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

  if(fl == NULL){
    if(mynode() == 0)printf("%s: Called with a null pointer\n", myname);
    terminate(1);
  }

  if( phases_in != 1){
    if(mynode() == 0)printf("BOTCH: %s needs phases in\n",myname); 
    terminate(1);
  }
  
  restore_qop_hisq_links_t(fl->flg, precision, links, &fl->options);
}

/*----------------------------------------*/
/* Accessors                              */
/*----------------------------------------*/

/* Return a list of fn_links types, one for each Naik epsilon */
fn_links_qop_t **
get_fm_links(fermion_links_t *fl){
  if(fl == NULL)return NULL;
  return get_qop_hisq_links_fn(fl->flg);
}

/* Return an fn_links type */
fn_links_qop_t *
get_fn_deps_links(fermion_links_t *fl){
  if(fl == NULL)return NULL;
  return get_qop_hisq_links_fn_deps(fl->flg);
}

/* Return the full FermionLinksHisq structure */
QOP_F3_FermionLinksHisq *
get_F_hisq_links(fermion_links_t *fl){
  if(fl == NULL)return NULL;
  return get_F_qop_hisq_links_hl(fl->flg);
}

/* Return the full FermionLinksHisq structure */
QOP_D3_FermionLinksHisq *
get_D_hisq_links(fermion_links_t *fl){
  if(fl == NULL)return NULL;
  return get_D_qop_hisq_links_hl(fl->flg);
}

QOP_hisq_coeffs_t *
get_action_coeffs_hisq(fermion_links_t *fl){
  if(fl == NULL)return NULL;
  return get_qop_hisq_links_ac(fl->flg);
}

int
get_n_naiks_hisq(fermion_links_t *fl){
  QOP_hisq_coeffs_t * ac;
  if(fl == NULL)return -1;
  ac = get_action_coeffs_hisq(fl);
  return get_n_naiks_qop(ac);
}

double *
get_eps_naik_hisq(fermion_links_t *fl){
  QOP_hisq_coeffs_t * ac;
  if(fl == NULL)return NULL;
  ac = get_action_coeffs_hisq(fl);
  return get_eps_naik_qop(ac);
}

int 
valid_fermion_links(fermion_links_t *fl, int precision){
  if(fl == NULL)return 0;
  return valid_qop_hisq_links(fl->flg, precision);
}

char *
get_action_parameter_string(fermion_links_t *fl){
  QOP_hisq_coeffs_t *ac = get_action_coeffs_hisq(fl);
  char *str = get_ac_string_hisq_qop(ac);

  return str;
}

