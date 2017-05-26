#include <iostream>
#include <cstdio>

#include <qphix/qphix_config.h>
#include <qphix/print_utils.h>
#include <qphix/geometry.h>

#include "../include/qphixj/qphixj.h"
#include "../include/qphixj/qphixj_arch.h"
#include <omp.h>

extern "C" {
#include "generic_includes.h"
#include "../include/generic_qphixj.h"
}

QPHIXJ_vec_layout_t vecLayout;

using namespace std;
using namespace QPhiX;

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

void QPHIXJ_init(QPHIXJ_vec_layout_t *user_vecLayout){

  // Copy vector layout
  vecLayout = *user_vecLayout;
  const int *machsize = get_logical_dimensions();

  // Work out local lattice size
  int subLattSize[4] = { nx/machsize[0], ny/machsize[1], nz/machsize[2], nt/machsize[3] };

  // Size of checkerboarded X-dimension

  int lX1h = subLattSize[0]/2;
  int lY = subLattSize[1];
  int lZ = subLattSize[2];
  int lT = subLattSize[3];
  
  int S = QPHIX_SOALEN;
  int V = VECLEN_DP;

  int By = vecLayout.By;
  int Bz = vecLayout.Bz;
  int PadXY = vecLayout.PadXY;
  int PadXYZ = vecLayout.PadXYZ;
  int NCores = vecLayout.NCores;
  int Sy = vecLayout.Sy;
  int Sz = vecLayout.Sz;

  // Diagnostic information:

  masterPrintf("Declared topology is: %d %d %d %d\n", 
	   machsize[0], machsize[1], machsize[2], machsize[3]);

  masterPrintf("VECLEN_DP=%d VECLEN_SP=%d VECLEN_HP=%d SOALEN=%d\n", VECLEN_DP, VECLEN_SP, VECLEN_HP, S);
  masterPrintf("Global Lattice Size = ");
  masterPrintf(" %d\n %d\n %d\n %d\n", nx, ny, nz, nt);
  masterPrintf("\n");

  masterPrintf("Local Lattice Size = ");
  for(int mu=0; mu < 4; mu++){ 
    masterPrintf(" %d", subLattSize[mu]);
  }
  masterPrintf("\n");

  if( lX1h % S != 0){
    masterPrintf("QPhiX ERROR: X length %d after checkerboarding must be divisible by SOALEN %d\n",
		 lX1h, S);
    abort();
  }

  masterPrintf("Block Sizes: By= %d Bz=%d\n", By, Bz);
  masterPrintf("Cores = %d\n", NCores);
  masterPrintf("SMT Grid: Sy=%d Sz=%d\n", Sy, Sz);
  masterPrintf("Pad Factors: PadXY=%d PadXYZ=%d\n", PadXY, PadXYZ);

  masterPrintf("\nQPHIXJ is Initialized\n");

}

void QPHIXJ_finalize(void){
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

  QPHIXJ_status_t status = QPHIXJ_SUCCESS;

  static QPHIXJ_vec_layout_t vec_layout = QPHIXJ_VEC_LAYOUT_DEFAULT;
  
  if(is_qphixj_env_setup > 0){
    /* Finalize so we can then reinitialize */
    QPHIXJ_finalize();
  }

  vec_layout.NCores = numCores;
  vec_layout.MinCt = minCt;

  node0_printf("Initializing QPhiX JLab\n");
  node0_printf("NumCores = %d, ThreadsPerCore = %d, minCt = %d\n", numCores, threads_per_core, minCt);

  QPHIXJ_init(&vec_layout);

  if(status != QPHIXJ_SUCCESS){
    node0_printf("Error initializing Qphixj\n");
    fflush(stdout);
    terminate(1);
  }

  fflush(stdout);

  is_qphixj_env_setup = 1;
  return status;
}
