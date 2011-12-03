/******************** fermion_links_from_site.c ****************************/
/* MIMD version 7 */

/* Teporary routines until we have removed the gauge file from the
   site structure */

#include "generic_ks_includes.h"
#include "../include/fermion_links.h"

fermion_links_t *create_fermion_links_from_site(int prec, int n_naiks, double *eps_naik){
  su3_matrix *links;
  fermion_links_t *fl;

  links = create_G_from_site();

#if FERM_ACTION == HISQ
  fl = create_fermion_links_hisq(prec, n_naiks, eps_naik, phases_in, links);
#else
  fl = create_fermion_links(prec, phases_in, links);
#endif

  free(links);
  return fl;
}

void restore_fermion_links_from_site(fermion_links_t *fl, int prec){
  su3_matrix *links;

  if(valid_fermion_links(fl, prec))return;

  links = create_G_from_site();

#if FERM_ACTION == HISQ
  restore_fermion_links_hisq(fl, prec, phases_in, links);
#else
  restore_fermion_links(fl, prec, phases_in, links);
#endif

  free(links);
}
