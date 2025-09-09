/******************* fermion_force_hisq_multi_grid.c ************************/
/* For the Grid interface */
/* MIMD version 7 */

/* PLACEHOLDER */

#include "generic_ks_includes.h"
#include "../include/generic_grid.h"
#include "../include/generic_ks_grid.h"
#include "../include/mGrid/mGrid.h"
#include "../include/openmp_defs.h"
#include <string.h>

#define GRID_MAX_NAIK 4

/*--------------------------------------------------------------------*/

extern GRID_4Dgrid *grid_full;

static GRID_hisq_coeffs_t hc = GRID_HISQ_COEFFS_DEFAULT;  /* Unsupported at present */

/*--------------------------------------------------------------------*/

void 
fermion_force_multi_hisq_grid(info_t* info, int prec, Real eps, Real *residues, 
			      su3_vector **multi_x, int nterms,
			      fermion_links_t *fl)
{
  char myname[] = "fermion_force_multi_hisq_grid";
  GRID_info_t grid_info;

  if(! grid_initialized()){
    node0_printf("%s: FATAL Grid has not been initialized\n", myname);
    terminate(1);
  }

  if(prec != 2){
    node0_printf("%s: WARNING, Overriding precision request %d and using 2 (double).\n",
		 myname, prec);
  }
  
  /* Compute weights, following the QOP convention for factors of 2 */
  double *epsv = (double *)malloc(sizeof(double)*nterms);
  for(int i = 0; i < nterms; i++) epsv[i] = eps*residues[i];

  su3_matrix* momentum = (su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));

  /* Support only double-precision calculation */
  GRID_D3_hisq_force(&grid_info, fl, epsv, multi_x, n_orders_naik, momentum,
		     grid_full);
  // append result
  site *s; int i;
  FORALLSITES_OMP(i,s,){
    anti_hermitmat ah3;
    for(int dir=0; dir<4; ++dir){
      make_anti_hermitian(&momentum[4*i + dir], &ah3);
      s->mom[dir].m00im    += ah3.m00im;
      s->mom[dir].m11im    += ah3.m11im;
      s->mom[dir].m22im    += ah3.m22im;   
      s->mom[dir].m01.real += ah3.m01.real;
      s->mom[dir].m02.real += ah3.m02.real;
      s->mom[dir].m12.real += ah3.m12.real;
      s->mom[dir].m01.imag += ah3.m01.imag;
      s->mom[dir].m02.imag += ah3.m02.imag;
      s->mom[dir].m12.imag += ah3.m12.imag;
    }
  } END_LOOP_OMP

  free(epsv);
  return;
  
}

/* fermion_links_hisq_load_grid.c */
