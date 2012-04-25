/**************** hisq_links_qop.c *****************************/
/* MILC Version 7 */

/* Methods for the hisq_links_qop_t "class"  */

#include "generic_ks_includes.h"
#include "../include/hisq_links_qop.h"
#include "../include/fn_links.h"
#include "../include/ks_action_coeffs_qop.h"
#include "../include/generic_qop.h"
#include "../include/generic_qopqdp.h"

#ifdef FLTIME
static const char *qop_prec[2] = {"F", "D"};
#endif

hisq_links_qop_t *
create_hisq_links_qop(QOP_hisq_coeffs_t *hc, int precision,
		      su3_matrix *links, int want_deps, int want_aux){

  hisq_links_qop_t *hl;
  QOP_info_t info = {0., 0., 0, 0, 0};
  fsu3_matrix **raw_F;
  dsu3_matrix **raw_D;
  QOP_F3_GaugeField *g_F;
  QOP_D3_GaugeField *g_D;
  /* Initialize to default values */
  QOP_opt_t qop_hl_opt[2] = {
    {.tag = "want_deps",.value=0},
    {.tag = "want_aux",.value=1}
  };

  char myname[] = "create_hisq_links_qop";
  
  hl = (hisq_links_qop_t *)malloc(sizeof(hisq_links_qop_t));
  if(hl == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }
  
  hl->phase = create_link_phase_info();
  hl->hl_F = NULL;
  hl->hl_D = NULL;

  /* Set links options, overriding defaults */
  qop_hl_opt[0].value = want_deps;
  qop_hl_opt[1].value = want_aux;

  QOP_hisq_links_set_opts(qop_hl_opt, 2);

  if(precision == 1){

    raw_F = create_raw4_F_G_from_field(links, EVENANDODD);
    g_F = QOP_F3_create_G_from_raw((float **)raw_F, QOP_EVENODD);
    destroy_raw4_F_G(raw_F);
    hl->hl_F = QOP_F3_hisq_create_L_from_G(&info, hc, g_F);
    QOP_F3_destroy_G(g_F);

  } else {

    raw_D = create_raw4_D_G_from_field(links, EVENANDODD);
    g_D = QOP_D3_create_G_from_raw((double **)raw_D, QOP_EVENODD);
    destroy_raw4_D_G(raw_D);
    hl->hl_D = QOP_D3_hisq_create_L_from_G(&info, hc, g_D);
    QOP_D3_destroy_G(g_D);

  }

#ifdef FLTIME
  node0_printf("FLTIME: time = %e (HISQ qop %s) flops/site = %d mflops = %e\n",
	       info.final_sec,qop_prec[precision-1],
	       (int)(info.final_flop*numnodes()/volume),
	       (Real)info.final_flop/(1e6*info.final_sec) );
#endif

  return hl;
}

/*-------------------------------------------------------------------*/
void 
destroy_hisq_links_qop(hisq_links_qop_t *hl){
  if(hl == NULL)return;

  destroy_link_phase_info(hl->phase);
  if(hl->hl_F != NULL)QOP_F3_hisq_destroy_L(hl->hl_F);
  if(hl->hl_D != NULL)QOP_D3_hisq_destroy_L(hl->hl_D);
  free(hl);
}

/*-------------------------------------------------------------------*/
/* Accessors                                                         */
/*-------------------------------------------------------------------*/

QOP_F3_FermionLinksHisq *
get_F_hisq_links_qop(hisq_links_qop_t *hl){
  return hl->hl_F;
}

QOP_D3_FermionLinksHisq *
get_D_hisq_links_qop(hisq_links_qop_t *hl){
  return hl->hl_D;
}

/*-------------------------------------------------------------------*/

/* Copy pointers and initialize phase info */

void
set_asqtad_links_from_hisq(fn_links_qop_t *fn, hisq_links_qop_t *hl, int i){
  //  char myname[] = "set_asqtad_links_from_hisq";

  if(fn == NULL)return;

  if(hl->hl_F != NULL)
    fn->al_F = QOP_F3_get_asqtad_links_from_hisq(hl->hl_F)[i];

  if(hl->hl_D != NULL)
    fn->al_D = QOP_D3_get_asqtad_links_from_hisq(hl->hl_D)[i];

  /* The global hl phase information has to be copied into fn, because
     there are multiple, independent fn links. */

  fn->phase = create_link_phase_info();
  copy_link_phase_info(fn->phase, hl->phase);
}

/* Clear copied pointers and destroy phase info */

void
unset_asqtad_links_from_hisq(fn_links_qop_t *fn){

  if(fn == NULL)return;

  free_fn_links_qop(fn);
  fn->al_D = NULL;
  fn->al_F = NULL;
  destroy_link_phase_info(fn->phase);
  fn->phase = NULL;
}

void
set_asqtad_deps_links_from_hisq(fn_links_qop_t *fn, hisq_links_qop_t *hl){

  if(fn == NULL)return;

  if(hl->hl_F != NULL)
    fn->al_F = QOP_F3_get_asqtad_deps_links_from_hisq(hl->hl_F);
  else
    fn->al_F = NULL;

  if(hl->hl_D != NULL)
    fn->al_D = QOP_D3_get_asqtad_deps_links_from_hisq(hl->hl_D);
  else
    fn->al_D = NULL;

  fn->phase = create_link_phase_info();
  copy_link_phase_info(fn->phase, hl->phase);
}


void
unset_asqtad_deps_links_from_hisq(fn_links_qop_t *fn){
  unset_asqtad_links_from_hisq(fn);
}

