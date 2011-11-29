/**************** fermion_links_eo_load_milc.c *****************************/
/* MILC Version 7 */

/* Methods for the trivial eo_links_t structure  */
/* Just copy the path table pointer */

#include "generic_ks_includes.h"
#include "../include/info.h"

/*-------------------------------------------------------------------*/
/* Fill in EO links                                           */
/*-------------------------------------------------------------------*/

/* This procedure is meant to be called only by fermion_links*.c
   routines. The public API is create_fermion_links and
   restore_fermion_links.  */

void
load_eo_links(info_t *info, eo_links_t *eo, ks_action_paths *ap, 
	      su3_matrix *links, int want_back){
  eo->ap = ap;
}
