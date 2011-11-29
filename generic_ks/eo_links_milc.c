/**************** eo_links_milc.c *****************************/
/* MILC Version 7 */
/* Methods for the FN links "class"  */

#include "generic_ks_includes.h"
#include "../include/eo_links.h"

/*-------------------------------------------------------------------*/
/* Create/destroy eo links                                           */
/*-------------------------------------------------------------------*/

/* The EO links are trivial.  They consist only of a copy of the path
   table pointer.  We follow the FN protool, so, here, we create only
   the structure and don't give a value to any member. */

eo_links_t *
create_eo_links(void){

  eo_links_t *eo;
  char myname[] = "create_eo_links";
  
  eo = (eo_links_t *)malloc(sizeof(eo_links_t));
  if(eo == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }
  
  eo->ap = NULL;
  return eo;
}

/*-------------------------------------------------------------------*/
void 
destroy_eo_links(eo_links_t *eo){
  if(eo == NULL)return;
  free(eo);
}

/*-------------------------------------------------------------------*/
/* Accessors                                                         */
/*-------------------------------------------------------------------*/

ks_action_paths *get_ap(eo_links_t *eo){
  return eo->ap;
}

