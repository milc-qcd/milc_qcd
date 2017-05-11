/******************* milc_to_qphixj_utilities.c ************************/
/* For the Qphixj interface */
/* MIMD version 7 */

/* 12/7/2016 Created by CD */

#include "../include/generic_qphixj.h"
#include "../include/generic.h"
#include <lattice.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>

static int is_qphixj_env_setup = 0;

/* The MILC layout routines list the coordinates and assume 4D */

static int milc_node_number(const int coords[]){
  return node_number(coords[0],coords[1],coords[2],coords[3]);
}

static int milc_node_index(const int coords[]){
  return node_index(coords[0],coords[1],coords[2],coords[3]);
}

QPHIXJ_evenodd_t milc2qphixj_parity(int milc_parity){
  switch(milc_parity){
  case(EVEN):       return QPHIXJ_EVEN;
  case(ODD ):       return QPHIXJ_ODD;
  case(EVENANDODD): return QPHIXJ_EVENODD;
  default:
    printf("milc2qphixj_parity: Bad MILC parity %d\n", milc_parity);
    terminate(1);
  }
  return (QPHIXJ_evenodd_t)-999;
}

static int qphixj2milc_parity(QPHIXJ_evenodd_t qphixj_parity){
  switch(qphixj_parity){
  case(QPHIXJ_EVEN):     return EVEN;
  case(QPHIXJ_ODD ):     return ODD;
  case(QPHIXJ_EVENODD ): return EVENANDODD;
  default:
    printf("qphixj2milc_parity: Bad QPHIXJ parity %d\n", qphixj_parity);
    terminate(1);
  }
  return -999;
}

QPHIXJ_status_t 
initialize_qphixj(void){
  /* \fixme - pass the qphixj params via setup */
  /* create mbench */

  int minCt = 1;
  int threads_per_core = 1;
  int numThreads = omp_get_max_threads();
  if(getenv("MINCT")) minCt = atoi(getenv("MINCT"));
  if(getenv("THREADS_PER_CORE")) threads_per_core = atoi(getenv("THREADS_PER_CORE"));
  int numCores = numThreads / threads_per_core;

  int status = 0;

  static int latsize[4];
  static QPHIXJ_layout_t layout;
  static QPHIXJ_vec_layout_t vec_layout = QPHIXJ_VEC_LAYOUT_DEFAULT;
  
  latsize[0] = nx;
  latsize[1] = ny;
  latsize[2] = nz;
  latsize[3] = nt;

  if(is_qphixj_env_setup > 0){
    /* Finalize so we can then reinitialize */
    finalize_qphixj();
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

  vec_layout.NCores = numCores;
  vec_layout.MinCt = minCt;

  node0_printf("Initializing Qphixj\n");
  node0_printf("NumCores = %d, ThreadsPerCore = %d, minCt = %d\n", numCores, threads_per_core, minCt);

  //status = QPHIXJ_init(&layout, precision, &vec_layout);
  QPHIXJ_init(&layout, &vec_layout);

  if(status){
    node0_printf("Error initializing Qphixj\n");
    fflush(stdout);
    terminate(1);
  }

  fflush(stdout);

  is_qphixj_env_setup = 1;
  return status;
}

void
finalize_qphixj(void){
  QPHIXJ_finalize();
}

