#ifndef _GENERIC_QUDA_H
#define _GENERIC_QUDA_H
/******************** generic_quda.h *********************************
*  MIMD version 7 	 				            *
*/

#ifdef HAVE_QUDA

#include "../include/openmp_defs.h"
#include <quda_milc_interface.h>

int initialize_quda(void);

static QudaMILCSiteArg_t newQudaMILCSiteArg() {
  QudaMILCSiteArg_t arg;
  arg.site = lattice;
  arg.link = NULL;
  arg.link_offset = (char*)lattice->link-(char*)lattice;
  arg.mom = NULL;
#ifdef MOM_SITE
  arg.mom_offset = (char*)lattice->mom-(char*)lattice;
#else
  arg.mom_offset = 0;
#endif
  arg.size = sizeof(*lattice);
  return arg;
}

void finalize_quda(void);

#include <string.h>

static inline void fast_copy(void *dest, const void *src, size_t n) {
  memcpy(dest, src, n);
}

// only use managed memory if enabled by QUDA
#ifndef USE_QUDA_MANAGED
#define qudaAllocateManaged qudaAllocatePinned
#define qudaFreeManaged qudaFreePinned
#endif

/*
  Allocate a pinned gauge-field array suitable for DMA transfer to the GPU
 */
static su3_matrix* create_G_quda(void) {
  return (su3_matrix*)qudaAllocatePinned(sites_on_node*4*sizeof(su3_matrix));
}

/*
  Extract the gauge field elements into a pinned array suitable for DMA transfer to the GPU
 */
static su3_matrix* create_G_from_site_quda(void) {
  su3_matrix *links = create_G_quda();
  int i;
  site *s;

  FORALLSITES_OMP(i,s,){
    fast_copy(links+4*i, s->link, 4*sizeof(su3_matrix));
  } END_LOOP_OMP

  return links;
}

/*
  Copy the momentum field elements into the site struct array
 */
static void copy_to_site_from_G_quda(su3_matrix *links) {
  int i;
  site *s;

  FORALLSITES_OMP(i,s,){
    fast_copy(s->link, links+4*i, 4*sizeof(su3_matrix));
  } END_LOOP_OMP
}

/*
  Free the pinned gauge-field array
 */
static void destroy_G_quda(su3_matrix *links) {
  qudaFreePinned(links);
}

/*
  Allocate a pinned momentum-field array suitable for DMA transfer to the GPU
 */
static anti_hermitmat* create_M_quda(void) {
  return (anti_hermitmat*)qudaAllocatePinned(sites_on_node*4*sizeof(anti_hermitmat));
}

#ifdef MOM_SITE
/*
  Extract the momentum field elements into a pinned array suitable for DMA transfer to the GPU
 */
static anti_hermitmat* create_M_from_site_quda(void) {
  anti_hermitmat* momentum = create_M_quda();
  int i;
  site *s;

  FORALLSITES_OMP(i,s,){
    fast_copy(momentum+4*i, s->mom, 4*sizeof(anti_hermitmat));
  } END_LOOP_OMP

  return momentum;
}

/*
  Copy the momentum field elements into the site struct array
 */
static void copy_to_site_from_M_quda(anti_hermitmat *momentum) {
  int i;
  site *s;

  FORALLSITES_OMP(i,s,){
    fast_copy(s->mom, momentum+4*i, 4*sizeof(anti_hermitmat));
  } END_LOOP_OMP
}
#endif

/*
  Free the pinned gauge-field array
 */
static void destroy_M_quda(anti_hermitmat *momentum) {
  qudaFreePinned(momentum);
}

#endif // HAVE_QUDA

#endif /* GENERIC_QUDA_H */
