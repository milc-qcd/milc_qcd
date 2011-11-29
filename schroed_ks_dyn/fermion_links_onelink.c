/******** fermion_links_onelink.c *********/

/* MIMD version 7 */

/* Dummy versions of standard generic_ks routines for the one-link KS
   action used in this application.  These procedures must be replaced
   if we upgrade to Asqtad.  They work only with the one-link
   d_congrad5.c and dslash.c, which do not use the fermion links
   structure. */

#include "schroed_ks_includes.h"

void load_longlinks(ferm_links_t *fn) {
  fn->lng = NULL;
}

void load_fatlinks(ferm_links_t *fn) {
  fn->fat = NULL;
}

void load_longbacklinks(ferm_links_t *fn){
  fn->lngback = NULL;
}

void load_fatbacklinks(ferm_links_t *fn){
  fn->fatback = NULL;
}

void load_ferm_links(ferm_links_t *fn){
  if(fn->valid == 1)return;

  load_fatlinks(fn);
  load_longlinks(fn);

#ifdef DBLSTORE_FN
  load_fatbacklinks(fn);
  load_longbacklinks(fn);
#endif

  fn->valid = 1;
}

