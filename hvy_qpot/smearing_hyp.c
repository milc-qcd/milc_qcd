/******************* smearing.c ****************************************/

/* MIMD version 7 */

/* Perform one iteration of "HYP 3D smearing" on the space-like links.
   Uses generic HYP-smearing routines with a flag that excludes TUP direction.
   A. Bazavov Oct 2012 */

#include "hvy_qpot_includes.h"
//#define TOL 1e-5
//#define MAXCOUNT 50

#include "../include/hyp_coeff.h"

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

// create a field with "n" entities on a site
// (n=4 is a usual link field)
static su3_matrix *create_mn_special(int n){
  char myname[] = "create_mn_special";
  su3_matrix *m;

  m = (su3_matrix *)special_alloc(sites_on_node*n*sizeof(su3_matrix));

  if(m==NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  return m;
}

static void destroy_mn_special(su3_matrix *m){
  special_free(m);
}


void smearing( void ) {

  su3_matrix *U_link, *hyp_link;
  int i,dir;
  site *st;
  hyp_coeffs_t hc;

  // create temporary storage
  U_link = create_mn_special(4);
  hyp_link = create_mn_special(4);

  // map original links into field
  FORALLSITES(i,st) {
    for(dir=XUP;dir<=TUP;dir++) {
      U_link[4*i+dir] = st->link[dir];
    }
  }

  // set parameters
  set_hyp_coeff( &hc, 0, hyp_alpha2, hyp_alpha3 );
  set_hyp_proj_method( &hc, HYP_SU3_TR_MAX, 9 );

  // HYP-smeared links, TUP excluded
  load_hyp_links(U_link, hyp_link, TUP, &hc);

  /* replace the original links, TUP links untouched */
  FORALLSITES(i,st)for(dir=XUP;dir<TUP;dir++){
    st->link[dir] = hyp_link[4*i+dir];
  }

  // free temporary storage
  destroy_mn_special(U_link);
  destroy_mn_special(hyp_link);

} /* smearing */
