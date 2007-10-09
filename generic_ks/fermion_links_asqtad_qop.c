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
load_fn_links(fn_links_t *fn, ks_action_paths *ap){

  if(PRECISION == 1)
    load_fn_links_F(fn, ap);
  else
    load_fn_links_D(fn, ap);
}

#ifdef DM_DU0
/* Wrappers for MILC call to QOP */
void load_fn_links_dmdu0(fn_links_t *fn, ks_action_paths *ap){

  if(PRECISION == 1)
    load_fn_links_dmdu0_F(fn, ap);
  else
    load_fn_links_dmdu0_D(fn, ap);
}
#endif

void
invalidate_fn_links(fn_links_t *fn)
{
  /* We must invalidate for both precisions */
  invalidate_fn_links_F(fn);
  invalidate_fn_links_D(fn);
}

void init_fn_links(fn_links_t *fn){
  fn->valid = 0;
  fn->fat = NULL;
  fn->lng = NULL;
  fn->fatback = NULL;
  fn->lngback = NULL;
}
