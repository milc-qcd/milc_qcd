/*********************** gauge_force_imp_grid.c  -- ****************************/
/* MIMD version 7 */
/* gauge action stuff for improved action
* T.D. and A.H. general gauge action updating code
* D.T. modified  5/97
* D.T. modified 12/97, optimized gauge_force a little
* D.T. modified 3/99, gauge action in include file
* C.D. split from gauge_stuff.c 10/06 */

/**#define GFTIME**/ /* For timing gauge force calculation */
#include "generic_includes.h"	/* definitions files and prototypes */
#include "../include/generic_grid.h"
#include "../include/openmp_defs.h"

extern GRID_4Dgrid *grid_full;

void imp_gauge_force_grid(Real eps, field_offset mom_off){
  char myname[] = "imp_gauge_force_grid";

  if(! grid_initialized()){
    node0_printf("%s: FATAL Grid has not been initialized\n", myname);
    terminate(1);
  }

#ifdef GFTIME
  int nflop = 153004;  /* For Symanzik1 action */
  double dtime = -dclock();
#endif

#if 0
  Real **loop_coeff = get_loop_coeff();
  const int num_loop_types = get_nloop();
  double *grid_loop_coeff = (double*)malloc(num_loop_types * sizeof(double));
#endif  

  site *s; int i;
  const Real eb3 = eps*beta/3.0;

  GRID_info_t grid_info;
  su3_matrix* momentum = (su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));
  su3_matrix* U = (su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));

  FORALLSITES_OMP(i,s,){
    int dir;
    FORALLUPDIR(dir){
      su3mat_copy(&s->link[dir], &U[4*i+dir]);
    }
  } END_LOOP_OMP;

  Real **loop_coeff = get_loop_coeff();
  int nloop = get_nloop();
  
  double cp  = loop_coeff[0][0];
  double cr  = (nloop > 1) ? -loop_coeff[1][0] : 0;
  double cpg = (nloop > 2) ? loop_coeff[2][0] : 0;

  if(MILC_PRECISION == 1)
    GRID_F3_gauge_force(&grid_info, U, eb3, cp, cr, cpg, momentum, grid_full);
  else
    GRID_D3_gauge_force(&grid_info, U, eb3, cp, cr, cpg, momentum, grid_full);

  // append result
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
  } END_LOOP_OMP;

  free(U);
  free(momentum);

#ifdef GFTIME
  dtime+=dclock();
  node0_printf("GFTIME:   time = %e (Symanzik1_GRID) mflops = %e\n",dtime,
	       nflop*(double)volume/(1e6*dtime*numnodes()) );
#endif

  return;
}

