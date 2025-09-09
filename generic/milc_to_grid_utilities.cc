/******************* milc_to_grid_utilities.cc ************************/
/* For the Grid interface */
/* MIMD version 7 */

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <Grid/Grid.h>
#include <Grid/communicator/Communicator.h>
#include <vector>
#include <iostream>
//#include <qmp.h>

#undef GRID_EXTERN

#include "../include/mGrid/mGrid_assert.h"

extern "C" {
  extern	int nx,ny,nz,nt;	/* lattice dimensions defined in lattice.h */
#include "../include/generic_grid.h"
  void check_create(void);
}
using namespace std;
using namespace Grid;
//using namespace Grid::QCD;

GRID_4Dgrid *grid_full;
GRID_4DRBgrid *grid_rb;

Coordinate squaresize;

static int grid_is_initialized = 0;

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

void finalize_grid()
{
  /* We omit MPI_Finalize() because most likely it will break a lot of things */
  Grid_unquiesce_nodes();
  GRID_destroy_4DRBgrid(grid_rb);
  GRID_destroy_4Dgrid(grid_full);
}

int grid_initialized(void){
  return grid_is_initialized;
}

#define MAXARG 13
#define MAXARGSTR 64

GRID_status_t 
initialize_grid(void){
  
  if(grid_is_initialized)
    return GRID_SUCCESS;

  /* We simulate the command line parameters to initialize Grid, and
   * call the standard routine. This doesn't look very elegant, but
   * this way we let Grid handle the parameters. */
  
  const int *machsize = get_logical_dimensions(); // dimensions of the array of MPI ranks
  squaresize.resize(4);
  squaresize[0] = nx/machsize[0];
  squaresize[1] = ny/machsize[1];
  squaresize[2] = nz/machsize[2];
  squaresize[3] = nt/machsize[3];

  int mpiX = machsize[0];
  int mpiY = machsize[1];
  int mpiZ = machsize[2];
  int mpiT = machsize[3];

  int  argc = MAXARG;
  char **argv = (char **)malloc(MAXARG*sizeof(char *));

  char tag_grid[] = "--grid";
  char val_grid[MAXARGSTR];
  snprintf (val_grid, MAXARGSTR, "%d.%d.%d.%d\0", nx, ny, nz, nt);

  char tag_mpi[] = "--mpi";
  char val_mpi[MAXARGSTR];
  snprintf (val_mpi, MAXARGSTR, "%d.%d.%d.%d\0",  mpiX,   mpiY,   mpiZ,   mpiT);

  char val_shm[MAXARGSTR];
#ifdef GRID_SHMEM_MAX
  snprintf (val_shm, MAXARGSTR, "%d\0", GRID_SHMEM_MAX);
  char tag_shm[] = "--shm";
#else
  snprintf (val_shm, MAXARGSTR, "");
  char tag_shm[] = "";
#endif

  char val_shmmpi[MAXARGSTR];
#ifdef GRID_SHMEM_MPI
  snprintf (val_shmmpi, MAXARGSTR, "%d\0", GRID_SHMEM_MPI);
  char tag_shmmpi[] = "--shm-mpi";
#else
  snprintf (val_shmmpi, MAXARGSTR, "");
  char tag_shmmpi[] = "";
#endif

  char val_device_mem[MAXARGSTR];
#ifdef GRID_DEVICE_MEM_MAX
  snprintf (val_device_mem, MAXARGSTR, "%d\0", GRID_DEVICE_MEM_MAX);
  char tag_device_mem[] = "--device-mem";
#else
  snprintf (val_device_mem, MAXARGSTR, "");
  char tag_device_mem[] = "";
#endif

  char tag_at[] = "--accelerator-threads";
  char val_at[MAXARGSTR];
#ifdef GRID_ACCELERATOR_THREADS
  int at = GRID_ACCELERATOR_THREADS;
#else
  int at = 8;  /* Default */
#endif
  snprintf (val_at, MAXARGSTR, "%d\0", at);
  
#ifdef GRID_COMMS_OVERLAP
  char tag_co[] = "--comms-overlap";
#else
  char tag_co[] = "";
#endif
  
  argv[0] = tag_grid; argv[1] = val_grid;
  argv[2] = tag_mpi;  argv[3] = val_mpi;
  argv[4] = tag_shm;  argv[5] = val_shm;
  argv[6] = tag_shmmpi; argv[7] = val_shmmpi;
  argv[8] = tag_device_mem;  argv[9] = val_device_mem;
  argv[10] = tag_at;   argv[11] = val_at;
  argv[12] = tag_co;

  if(mynode()==0)printf("Calling Grid_init with %s %s %s %s %s %s %s %s %s %s %s %s %s\n",argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9],argv[10],argv[11],argv[12]);
  Grid_init(&argc, &argv);

  grid_full = GRID_create_grid();
  grid_rb = GRID_create_RBgrid(grid_full);

  Coordinate latt_size   = GridDefaultLatt();
  Coordinate mpi_layout  = GridDefaultMpi();

  if(mynode()==0){
    printf("milc_to_grid_utilities: Initialized Grid with args\n%s %s\n%s %s\n%s %s\n",
	       argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]);
    printf("latt_size: %d %d %d %d\n", latt_size[0], latt_size[1], latt_size[2], latt_size[3]);
    printf("mpi_layout: %d %d %d %d\n", mpi_layout[0], mpi_layout[1], mpi_layout[2], mpi_layout[3]);
    printf("Grid threads %d\n", Grid::GridThread::GetThreads());
  }
  fflush(stdout);

  free(argv);

  grid_is_initialized = 1;
  return GRID_SUCCESS;

}

static Grid::CartesianCommunicator *grid_cart = NULL;

// Set up communicator for quick access to the rank mapping

void setup_grid_communicator(int peGrid[]){
  if(! grid_is_initialized){
    printf("setup_grid_communicator: Grid must first be initialized\n");
    terminate(1);
  }

  Coordinate processors;
  for(int i=0;i<4;i++) processors.push_back(peGrid[i]);
  grid_cart = new Grid::CartesianCommunicator(processors);
}

// Return the coordinate as assigned by Grid
// Code from CPS via Chulwoo Jung

int *query_grid_node_mapping(void){

  if(! grid_is_initialized){
    if(mynode()==0)printf("query_grid_node_mapping: Grid must first be initialized\n");
    terminate(1);
  }

  if(grid_cart == NULL){
    if(mynode()==0)printf("query_grid_node_mapping: Must call setup_grid_communicator first\n");
    terminate(1);
  }

  static int pePos[4];  /* Position of this process in the grid.*/

  for(int i=0;i<4;i++){
    pePos[i] = grid_cart->_processor_coor[i];
  }
  return pePos;
}

// The MPI rank that has the processor defined by lexrank
// COULD USE THIS WHEN GRIDS_COMMS_MPI3 IS DEFINED
// int grid_lexicographic_to_worldrank(int lexrank){
//   if(! grid_is_initialized){
//     printf("grid_lexicographic_to_worldrank: Grid must first be initialized\n");
//     terminate(1);
//   }
// 
//   if(grid_cart == NULL){
//     printf("grid_lexicographic_to_worldrank: Must call setup_grid_communicator first\n");
//     terminate(1);
//   }
// 
//   int worldrank = grid_cart->LexicographicToWorldRank[lexrank];
//   return worldrank;
// }

int grid_rank_from_processor_coor(int x, int y, int z, int t){
  if(! grid_is_initialized){
    if(mynode()==0)printf("grid_lexicographic_to_worldrank: Grid must first be initialized\n");
    terminate(1);
  }

  if(grid_cart == NULL){
    if(mynode()==0)printf("grid_lexicographic_to_worldrank: Must call setup_grid_communicator first\n");
    terminate(1);
  }

  Coordinate coor(4);
  coor[0] = x; coor[1] = y; coor[2] = z; coor[3] = t;
  int worldrank = grid_cart->RankFromProcessorCoor(coor);
  return worldrank;
}

void grid_coor_from_processor_rank(int coords[], int worldrank){
  if(! grid_is_initialized){
    if(mynode()==0)printf("grid_coor_from_processor_rank: Grid must first be initialized\n");
    terminate(1);
  }

  if(grid_cart == NULL){
    if(mynode()==0)printf("grid_coor_from_processor_rank: Must call setup_grid_communicator first\n");
    terminate(1);
  }

  Coordinate coor(4);
  grid_cart->ProcessorCoorFromRank(worldrank, coor);

  for(int i = 0; i < 4; i++)
    coords[i] = coor[i];

}

