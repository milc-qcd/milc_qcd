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
  fn->al_F_allocated = 0;
  fn->al_F = NULL;
  fn->al_D = NULL;
  fn->fat = NULL;
  fn->lng = NULL;

  return fn;
}

void
load_fn_links_qop(fn_links_qop_t *fn, QOP_asqtad_coeffs_t *coeffs, 
		  int precision, su3_matrix *links, int want_back){
  
  QOP_info_t info = {0., 0., 0, 0, 0};
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
  if(fn->fat != NULL)free(fn->fat);
  if(fn->lng != NULL)free(fn->lng);
  free(fn);
}

/*-------------------------------------------------------------------*/

QOP_F3_FermionLinksAsqtad *
get_F_asqtad_links(fn_links_qop_t *fn){
  /* If we have only the double-precision links, copy them to single precision */
  if(fn->al_F == NULL && fn->al_D != NULL){
    /* Here we actually allocate space for the al_F links */
    fn->al_F = QOP_FD3_asqtad_create_L_from_L(fn->al_D);
    fn->al_F_allocated = 1;
  }
  return fn->al_F;
}

QOP_D3_FermionLinksAsqtad *
get_D_asqtad_links(fn_links_qop_t *fn){
  return fn->al_D;
}

void
free_fn_links_qop(fn_links_qop_t *fn){
  /* Undo the copy that we did in get_F_asqtad_links above */
  if(fn->al_F_allocated)
    QOP_F3_asqtad_destroy_L(fn->al_F);
  fn->al_F = NULL;
  fn->al_F_allocated = 0;
}

/*-------------------------------------------------------------------*/

/* Extract MILC-style fat and long links from the QOP FN members of
   the fn_links_qop_t structure. The extracted values are read-only.
   That is, we don't support changing the fat and long links and
   expect the changes to propagate back to the QOP FN links */

/* Policy: We might not have any long links here, so we don't create
   space for them or extract them in get_fatlinks -- only if
   explicitly requested through the get_lnglinks call.  But we always
   have fat links, so we create space for fat links and extract them
   in both get_fatlinks and get_lnglinks. */


su3_matrix *get_fatlinks(fn_links_qop_t *fn){
  
  QOP_D3_FermionLinksAsqtad *fn_D;
  QOP_F3_FermionLinksAsqtad *fn_F;

  /* If the fat links have already been extracted, use them */
  if(fn->fat != NULL)
    return fn->fat;

  if(fn->lng != NULL){
    node0_printf("get_lnglinks: unexpected fat==NULL lng!=NULL\n");
    terminate(1);
  }

  fn->fat = create_G(); /* 4 su3_matrices per site */

  /* Try the double-precision values first */
  fn_D = get_D_asqtad_links(fn);
  if(fn_D != NULL){
    /* Unloads only fat because lng == NULL */
    unload_D_L_to_fields( fn->fat, fn->lng, fn_D, EVENANDODD);
    return fn->fat;
  }

  fn_F = get_F_asqtad_links(fn);
  if(fn_F != NULL){
    /* Unloads only fat because lng == NULL */
    unload_F_L_to_fields( fn->fat, fn->lng, fn_F, EVENANDODD);
    return fn->fat;
  }
  
  return NULL;
}

su3_matrix *get_lnglinks(fn_links_qop_t *fn){
  
  QOP_D3_FermionLinksAsqtad *fn_D;
  QOP_F3_FermionLinksAsqtad *fn_F;

  /* If the lng links have already been extracted, use them */
  if(fn->lng != NULL)
    return fn->lng;

  if(fn->fat == NULL)
    fn->fat = create_G(); /* 4 su3_matrices per site */

  fn->lng = create_G();

  /* Try the double-precision values first */
  fn_D = get_D_asqtad_links(fn);
  if(fn_D != NULL){
    /* Unloads both fat and lng */
    unload_D_L_to_fields( fn->fat, fn->lng, fn_D, EVENANDODD);
    return fn->lng;
  }

  fn_F = get_F_asqtad_links(fn);
  if(fn_F != NULL){
    /* Unloads both fat and lng */
    unload_F_L_to_fields( fn->fat, fn->lng, fn_F, EVENANDODD);
    return fn->lng;
  }
  
  return NULL;
}

/* We don't support the backward link extractions for now.  The
   routines that need them work directly with the fn_links_qop_t
   structure instead. */

su3_matrix *get_fatbacklinks(fn_links_qop_t *fn){
  return NULL;
}

su3_matrix *get_lngbacklinks(fn_links_qop_t *fn){
  return NULL;
}




