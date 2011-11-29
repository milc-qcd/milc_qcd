/**************** fn_links_qop.c *****************************/
/* MILC Version 7 */

/* Methods for the fn_links_qop_t "class"  */

#include "generic_ks_includes.h"
#include "../include/fn_links_qop.h"
#include "../include/ks_action_coeffs_qop.h"
#include "../include/generic_qop.h"
#include "../include/generic_qopqdp.h"


fn_links_qop_t *
create_fn_links_qop(void){

  fn_links_qop_t *fn;
  char myname[] = "create_fn_links_qop";
  
  fn = (fn_links_qop_t *)malloc(sizeof(fn_links_qop_t));
  if(fn == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }
  
  fn->phase = NULL;
  fn->al_F = NULL;
  fn->al_D = NULL;

  return fn;
}

void
load_fn_links_qop(fn_links_qop_t *fn, QOP_asqtad_coeffs_t *coeffs, 
		  int precision, su3_matrix *links, int want_back){
  
  QOP_info_t info;
  fsu3_matrix **raw_F;
  dsu3_matrix **raw_D;
  QOP_F3_GaugeField *g_F;
  QOP_D3_GaugeField *g_D;
  //  char myname[] = "load_fn_links_qop";

  destroy_link_phase_info(fn->phase);
  fn->phase = create_link_phase_info();

  if(precision == 1){

    raw_F = create_raw4_F_G_from_field(links, EVENANDODD);
    g_F = QOP_F3_create_G_from_raw((float **)raw_F, QOP_EVENODD);
    destroy_raw4_F_G(raw_F);
    fn->al_F = QOP_F3_asqtad_create_L_from_G(&info, coeffs, g_F);
    QOP_F3_destroy_G(g_F);

  } else {

    raw_D = create_raw4_D_G_from_field(links, EVENANDODD);
    g_D = QOP_D3_create_G_from_raw((double **)raw_D, QOP_EVENODD);
    destroy_raw4_D_G(raw_D);
    fn->al_D = QOP_D3_asqtad_create_L_from_G(&info, coeffs, g_D);
    QOP_D3_destroy_G(g_D);

  }

}

/*-------------------------------------------------------------------*/
void 
destroy_fn_links_qop(fn_links_qop_t *fn){
  if(fn == NULL)return;

  destroy_link_phase_info(fn->phase);
  if(fn->al_F != NULL)QOP_F3_asqtad_destroy_L(fn->al_F);
  if(fn->al_D != NULL)QOP_D3_asqtad_destroy_L(fn->al_D);
  free(fn);
}

/*-------------------------------------------------------------------*/

QOP_F3_FermionLinksAsqtad *
get_F_asqtad_links(fn_links_qop_t *fn){
  return fn->al_F;
}

QOP_D3_FermionLinksAsqtad *
get_D_asqtad_links(fn_links_qop_t *fn){
  return fn->al_D;
}

