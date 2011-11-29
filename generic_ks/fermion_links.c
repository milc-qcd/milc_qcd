/****** fermion_links.c  -- ******************/
/* MIMD version 7 */

/* March 2011 CD */

/* Options 

   Design considerations for options:

   There are a few global options that affect the building of the
   links structure.  They are set by calls to the routines below
   before first calling to create the links.  When the links are
   created, the prevailing options are copied into the links structure
   where they are subsequently unchangeable.

 */

/* Defaults ...

   No du0
   No deps
   Yes we want the HISQ auxiliaries for a fermion force calculation.
   No backward links

*/

#include <stdlib.h>
#include <stdio.h>
#include "../include/fermion_links.h"
#ifdef HAVE_QOP
#include "../include/fermion_links_qop.h"
#else
#include "../include/fermion_links_milc.h"
#endif
#include "../include/comdefs.h"

static ferm_links_options_t options = {0, 0, 1, 0};

/*---------------------------------------------------------------*/
/* Set global option flags.  Must be set BEFORE creating the links
   structure */
/*---------------------------------------------------------------*/

void fermion_links_want_du0(int request){
  options.want_du0 = request;
}

void fermion_links_want_deps(int request){
  options.want_deps = request;
}

void fermion_links_want_aux(int request){
  options.want_aux = request;
}

void fermion_links_want_back(int request){
  options.want_back = request;
}

/*----------------------------------------*/
/* Copy the options values                */
/*----------------------------------------*/

static void load_fl_options(ferm_links_options_t *dest){
  dest->want_du0 = options.want_du0;
  dest->want_deps = options.want_deps;
  dest->want_aux = options.want_aux;
  dest->want_back = options.want_back;
}

/*----------------------------------------*/
/* Create the fermion links structure */
/* This procedure is "private".  Use
   create_fermion_links_asqtad or
   create_fermion_links_hisq as appropriate */
/*----------------------------------------*/

/* To be called only by fermion_links*.c routines */
fermion_links_t *create_fermion_links_t(void){

  fermion_links_t *fl;
  char myname[] = "create_fermion_links";

  fl = (fermion_links_t *)malloc(sizeof(fermion_links_t));
  if(fl == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  load_fl_options(&fl->options);

  fl->flg = NULL;

  return fl;
}

/*----------------------------------------*/
/* Destroy the fermion links structure */
/*----------------------------------------*/

/* To be called only by fermion_links*.c routines */
void destroy_fermion_links_t(fermion_links_t *fl){

  if(fl == NULL)return;
  free(fl);
}

int fermion_links_get_n_naiks(fermion_links_t *fl){
  if(fl == NULL)return 0;

#if FERM_ACTION == HISQ
  return get_n_naiks_hisq(fl);
#else
  return 1;
#endif

}

double *fermion_links_get_eps_naik(fermion_links_t *fl){
  if(fl == NULL)return NULL;

#if FERM_ACTION == HISQ
  return get_eps_naik_hisq(fl);
#else
  return NULL;
#endif

}
