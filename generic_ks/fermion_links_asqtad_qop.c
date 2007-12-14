/****************** fermion_links_asqtad_qop.c ***********************/
/* MIMD version 7 */

/* This is the MILC wrapper for SciDAC Level 3 QOP link smearing */
/* These are generic entry points taking the prevailing MILC precision */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/generic_qop.h"
#include "../include/generic_ks_qop.h"


/*********************************************************************/
/* Create fat and long links and qop_links                           */
/*********************************************************************/
/* Wrappers for MILC call to QOP */
void 
load_ferm_links(ferm_links_t *fn, ks_action_paths *ap){

  if(PRECISION == 1)
    load_ferm_links_F(fn, ap);
  else
    load_ferm_links_D(fn, ap);
}

#ifdef DM_DU0
/* Wrappers for MILC call to QOP */
void load_ferm_links_dmdu0(ferm_links_t *fn, ks_action_paths *ap){

  if(PRECISION == 1)
    load_ferm_links_dmdu0_F(fn, ap);
  else
    load_ferm_links_dmdu0_D(fn, ap);
}
#endif

void
invalidate_all_ferm_links(ferm_links_t *fn)
{
  /* We must invalidate for both precisions */
  invalidate_all_ferm_links_F(fn);
  invalidate_all_ferm_links_D(fn);
}

void init_ferm_links(ferm_links_t *fn){
  fn->valid = 0;
  fn->fat = NULL;
  fn->lng = NULL;
  fn->fatback = NULL;
  fn->lngback = NULL;
  fn->ap = NULL;
  fn->valid_qop_F = 0;
  fn->valid_qop_D = 0;
  fn->qop_F_l = NULL;
  fn->qop_D_l = NULL;
}
