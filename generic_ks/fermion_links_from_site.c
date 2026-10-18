/******************** fermion_links_from_site.c ****************************/
/* MIMD version 7 */

/* Temporary routines until we have removed the gauge file from the
   site structure */

#include "generic_ks_includes.h"
#include "../include/fermion_links.h"
#include "../include/openmp_defs.h"

#if defined(USE_FL_GPU) && defined(HAVE_QUDA)
#include "../include/generic_quda.h"
#endif

/*----------------------------------------------------------------------*/
/* Multiply a copy of the thin links by the U(1) phase
   exp( i * charge * A_mu(x) ).

   The links array is a private copy made from the site structure, so
   the thin links in lattice[i].link are left untouched.  A NULL u1_A
   or a zero charge is a no-op, which lets the ordinary entry points
   below share this code path. */

static void
mult_links_by_u1_phase(su3_matrix *links, Real charge, Real *u1_A){
  int i, dir;
  site *s;

  if(u1_A == NULL || charge == 0.)return;

  FORALLSITES_OMP(i,s,private(dir)){
    FORALLUPDIR(dir){
      complex z = ce_itheta( charge*u1_A[4*i+dir] );   /* Qe * A(mu,x) */
      su3_matrix tmp = links[4*i+dir];
      /* not done in place: CMUL cannot alias its operands */
      c_scalar_mult_su3mat( &tmp, &z, links + 4*i + dir );
    }
  } END_LOOP_OMP;
}

/*----------------------------------------------------------------------*/
/* Copy the thin links out of the site structure, apply the U(1) phase
   and build the fermion links from the product. */

fermion_links_t *
create_fermion_links_from_site_u1(int prec, int n_naiks, double *eps_naik,
				  Real charge, Real *u1_A){
  fermion_links_t *fl;

#if defined(USE_FL_GPU) && defined(HAVE_QUDA)
  su3_matrix *links = create_G_from_site_quda();
#else
  su3_matrix *links = create_G_from_site();
#endif

  mult_links_by_u1_phase(links, charge, u1_A);

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

/*----------------------------------------------------------------------*/
void
restore_fermion_links_from_site_u1(fermion_links_t *fl, int prec,
				   Real charge, Real *u1_A){

  if(valid_fermion_links(fl, prec))return;

#if defined(USE_FL_GPU) && defined(HAVE_QUDA)
  su3_matrix *links = create_G_from_site_quda();
#else
  su3_matrix *links = create_G_from_site();
#endif

  mult_links_by_u1_phase(links, charge, u1_A);

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

/*----------------------------------------------------------------------*/
/* The matching destructor, so a caller need not repeat the HISQ test */

void
destroy_fermion_links_from_site(fermion_links_t *fl){
#if FERM_ACTION == HISQ
  destroy_fermion_links_hisq(fl);
#else
  destroy_fermion_links(fl);
#endif
}

/*----------------------------------------------------------------------*/
/* The plain entry points: no U(1) field */

fermion_links_t *create_fermion_links_from_site(int prec, int n_naiks, double *eps_naik){
  return create_fermion_links_from_site_u1(prec, n_naiks, eps_naik, 0., NULL);
}

void restore_fermion_links_from_site(fermion_links_t *fl, int prec){
  restore_fermion_links_from_site_u1(fl, prec, 0., NULL);
}
