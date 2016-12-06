/******************* milc_to_grid_utilities.c ************************/
/* For the Grid interface */
/* MIMD version 7 */

/* 12/03/16 Created by CD */

#include "../include/generic_grid.h"
#include "../include/generic.h"
#include <lattice.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>

static int is_grid_env_setup = 0;

/* The MILC layout routines list the coordinates and assume 4D */

static int milc_node_number(const int coords[]){
  return node_number(coords[0],coords[1],coords[2],coords[3]);
}

static int milc_node_index(const int coords[]){
  return node_index(coords[0],coords[1],coords[2],coords[3]);
}

GRID_evenodd_t milc2grid_parity(int milc_parity){
  switch(milc_parity){
  case(EVEN):       return GRID_EVEN;
  case(ODD ):       return GRID_ODD;
  case(EVENANDODD): return GRID_EVENODD;
  default:
    printf("milc2grid_parity: Bad MILC parity %d\n", milc_parity);
    terminate(1);
  }
  return (GRID_evenodd_t)-999;
}

int grid2milc_parity(GRID_evenodd_t grid_parity){
  switch(grid_parity){
  case(GRID_EVEN):     return EVEN;
  case(GRID_ODD ):     return ODD;
  case(GRID_EVENODD ): return EVENANDODD;
  default:
    printf("grid2milc_parity: Bad GRID parity %d\n", grid_parity);
    terminate(1);
  }
  return -999;
}

GRID_status_t 
initialize_grid(int precision){
  /* \fixme - pass the grid params via setup */
  /* create mbench */

  int minCt = 1;
  int threads_per_core = 1;
  int numThreads = omp_get_max_threads();
  if(getenv("MINCT")) minCt = atoi(getenv("MINCT"));
  if(getenv("THREADS_PER_CORE")) threads_per_core = atoi(getenv("THREADS_PER_CORE"));
  int numCores = numThreads / threads_per_core;

  int status = 0;

  static int latsize[4];
  static GRID_layout_t layout;
  
  latsize[0] = nx;
  latsize[1] = ny;
  latsize[2] = nz;
  latsize[3] = nt;

  if(is_grid_env_setup > 0){
    if(is_grid_env_setup == precision)
      return status;
    else 
      /* Finalize so we can then initialize with the new precision */
      finalize_grid();
  }

  layout.node_number = milc_node_number;
  layout.node_index = milc_node_index;
  layout.latdim = 4;
  layout.latsize = latsize;
  layout.machdim = 4;
  layout.machsize = (int *)get_logical_dimensions();
  layout.this_node = this_node;
  layout.even_sites_on_node = even_sites_on_node;
  layout.sites_on_node = sites_on_node;

  node0_printf("Initializing Grid for precision %d\n", precision);
  node0_printf("NumCores = %d, ThreadsPerCore = %d, minCt = %d\n", numCores, threads_per_core, minCt);

  status = GRID_init(&layout, precision);

  if(status){
    node0_printf("Error initializing Grid\n");
    fflush(stdout);
    terminate(1);
  }

  fflush(stdout);

  is_grid_env_setup = precision;
  return status;
}

void
finalize_grid(void){
  GRID_finalize();
}


