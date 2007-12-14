/****** fermion_links_eo.c  -- ******************/
/* MIMD version 7 */

/* External entry points

   init_ferm_links
   load_ferm_links
   load_ferm_links_dmdu0 (ifdef DM_DU0)
   invalidate_all_ferm_links

 */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#define IMP_QUARK_ACTION_INFO_ONLY
#include <quark_action.h>

void load_ferm_links(ferm_links_t *fn, ks_action_paths *ap){

  if(fn->valid == 1)return;
  fn->ap = ap;
  fn->valid = 1;
}

#ifdef DM_DU0
void load_ferm_links_dmdu0(ferm_links_t *fn, ks_action_paths *ap){
  if(fn->valid == 1)return;
  fn->ap = ap;
  fn->valid = 1;
}
#endif

void
invalidate_all_ferm_links(ferm_links_t *fn)
{
  fn->valid = 0;
}

void init_ferm_links(ferm_links_t *fn){
  fn->valid = 0;
  fn->fat = NULL;
  fn->lng = NULL;
  fn->fatback = NULL;
  fn->lngback = NULL;
  fn->ap = NULL;
}
