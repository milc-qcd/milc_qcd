/******************** fermion_links_from_site.c ****************************/
/* MIMD version 7 */

/* Temporary routines until we have removed the gauge file from the
   site structure */

#include "generic_ks_includes.h"
#include "../include/fermion_links.h"

#if defined(USE_FL_GPU) && defined(HAVE_QUDA)
#include "../include/generic_quda.h"
#endif

fermion_links_t *create_fermion_links_from_site(int prec, int n_naiks, double *eps_naik){
  fermion_links_t *fl;

#if defined(USE_FL_GPU) && defined(HAVE_QUDA)
  su3_matrix *links = create_G_from_site_quda();
#else
  su3_matrix *links = create_G_from_site();
#endif

#if FERM_ACTION == HISQ
  fl = create_fermion_links_hisq(prec, n_naiks, eps_naik, phases_in, links);
#else
  fl = create_fermion_links(prec, phases_in, links);
#endif

#if defined(USE_FL_GPU) && defined(HAVE_QUDA)
  destroy_G_quda(links);
#else
  free(links);
#endif
  return fl;
}

void restore_fermion_links_from_site(fermion_links_t *fl, int prec){

  if(valid_fermion_links(fl, prec))return;

#if defined(USE_FL_GPU) && defined(HAVE_QUDA)
  su3_matrix *links = create_G_from_site_quda();
#else
  su3_matrix *links = create_G_from_site();
#endif

#if FERM_ACTION == HISQ
  restore_fermion_links_hisq(fl, prec, phases_in, links);
#else
  restore_fermion_links(fl, prec, phases_in, links);
#endif

#if defined(USE_FL_GPU) && defined(HAVE_QUDA)
  destroy_G_quda(links);
#else
  free(links);
#endif
}
