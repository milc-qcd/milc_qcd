/****** imp_gauge_action_grid.c  -- ******************/

/* MIMD version 7 */
/* gauge action stuff for improved action
* T.D. and A.H. general gauge action updating code
* D.T. modified  5/97
* D.T. modified 12/97, optimized gauge_force a little
* D.T. modified 3/99, gauge action in include file
* E.W. modified 7/22, split off gauge action, added GPU implementation */

#include "generic_includes.h"
#include "../include/generic_grid.h"
#include "../include/openmp_defs.h"
extern GRID_4Dgrid *grid_full;

double imp_gauge_action_grid() {
  char myname[] = "imp_gauge_action_grid";
  
  if(! grid_initialized()){
    node0_printf("%s: FATAL Grid has not been initialized\n", myname);
    terminate(1);
  }

#if 0  
  /* these are for loop_table  */
  int ln,iloop;
  
  /* get loop variables from functions */
  const int max_length = get_max_length();
  const int nloop = get_nloop();
  const int nreps = get_nreps();
  const int *loop_length = get_loop_length();
  const int *loop_num = get_loop_num();
  int ***loop_table = get_loop_table();
  Real **loop_coeff_milc = get_loop_coeff();
  
  if (nreps != 1){
    printf("imp_gauge_action_gpu: Does not support nreps != 1, disable gauge action offload\n");
    terminate(1);
  }

  // Count total number of loops
  int num_paths = 0;
  for (iloop = 0; iloop < nloop; iloop++)
    for (ln = 0; ln < loop_num[iloop]; ln++)
      num_paths++;
  
#ifdef GATIME
  int nlinks = 0;
  for (iloop = 0; iloop < nloop; iloop++) {
    nlinks += loop_num[iloop] * loop_length[iloop];
  }
  int nflop = 198 * nlinks + 8 * num_paths; /* For any action */
  double dtime = -dclock();
#endif
  
  // Storage for input paths
  int **input_path_buf = (int**)malloc(num_paths * sizeof(int*));
  for (i = 0; i < num_paths; i++)
    input_path_buf[i] = (int*)malloc(max_length * sizeof(int));
  
  // Storage for path lengths
  int *path_length = (int*)malloc(num_paths * sizeof(int));
  
  // Storage for loop coefficients
  double *loop_coeff = (double*)malloc(num_paths * sizeof(double));
  
  // Populate arrays
  num_paths = 0;
  for (iloop = 0; iloop < nloop; iloop++) {
    length = loop_length[iloop];
    for (ln = 0; ln < loop_num[iloop]; ln++) {
      path_length[num_paths] = length; // path length
      loop_coeff[num_paths] = 1.0; // due to the "3. - [...]" convention below, we'll wait to scale then
      for (i = 0; i < length; i++)
	input_path_buf[num_paths][i] = loop_table[iloop][ln][i];
      num_paths++;
    }
  }

#else
  int nflop = 0;
#endif

  int i; site *s;
  int total_dyn_flavors = 0;
  for(i = 0; i < n_dyn_masses; i++){
    total_dyn_flavors += dyn_flavors[i];
  }

  const Real eb3 = 1.;
  
  GRID_info_t grid_info;
  su3_matrix* momentum = (su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));
  su3_matrix* U = (su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));

  FORALLSITES_OMP(i,s,){
    int dir;
    FORALLUPDIR(dir){
      su3mat_copy(&s->link[dir], &U[4*i+dir]);
    }
  } END_LOOP_OMP;

  double g_action;
  if(MILC_PRECISION == 1){
    g_action = GRID_F3_gauge_action(&grid_info, U, eb3, u0, total_dyn_flavors, grid_full);
  } else {
    g_action = GRID_D3_gauge_action(&grid_info, U, eb3, u0, total_dyn_flavors, grid_full);
  }
  
#if 0
  free(loop_coeff);
  free(path_length);
  for (i = 0; i < num_paths; i++)
    free(input_path_buf[i]);
  free(input_path_buf);
#endif
  
  free(U);

#ifdef GATIME
  dtime+=dclock();
  node0_printf("GATIME:   time = %e (GaugeAction_Symanzik1_GRID) mflops = %e\n",dtime,
	       nflop*(double)volume/(1e6*dtime*numnodes()) );
#endif

  return( g_action );
} /* imp_gauge_action_grid */


