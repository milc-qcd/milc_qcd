/******** layout_qdp.c *********/
/* MIMD version 7 */
/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES */
/* This version uses QDP for layout */
/*
   setup_layout() does any initial setup.  When it is called the
     lattice dimensions nx,ny,nz and nt have been set.
     This routine sets the global variables "sites_on_node",
     "even_sites_on_node" and "odd_sites_on_node".
   node_number(x,y,z,t) returns the node number on which a site lives.
   node_index(x,y,z,t) returns the index of the site on the node - ie the
     site is lattice[node_index(x,y,z,t)].
   get_logical_dimensions() returns the machine dimensions
   get_logical_coordinates() returns the mesh coordinates of this node
   These routines will change as we change our minds about how to distribute
     sites among the nodes.  Hopefully the setup routines will work for any
     consistent choices. (ie node_index should return a different value for
     each site on the node.)
*/
#include "generic_includes.h"
#include <qdp.h>
#include <qmp.h>

static const int* dim_mach;

static int nodes_per_ionode[4];
static int *ionodegeomvals = NULL; /* ionode partitions */

static void lex_coords(int coords[], const int dim, const int size[], 
	   const size_t rank)
{
  int d;
  size_t r = rank;

  for(d = 0; d < dim; d++){
    coords[d] = r % size[d];
    r /= size[d];
  }
}

/*------------------------------------------------------------------*/
/* Convert coordinate to linear lexicographic rank (inverse of
   lex_coords) */

static size_t lex_rank(const int coords[], int dim, const int size[])
{
  int d;
  size_t rank = coords[dim-1];

  for(d = dim-2; d >= 0; d--){
    rank = rank * size[d] + coords[d];
  }
  return rank;
}

#ifdef FIX_IONODE_GEOM

/*------------------------------------------------------------------*/
/* Initialize io_node function */

static void init_io_node(){
  int i;
  int status = 0;

  if(ionodegeom() == NULL){
    ionodegeomvals = ionode_geometry;
    node0_printf("Setting ionodegeomvals to %d %d %d %d\n",
		 ionodegeomvals[0], ionodegeomvals[1],
		 ionodegeomvals[2], ionodegeomvals[3]);
  } else {
    node0_printf("init_io_node: Command line ionode geometry overrides request\n");
    ionodegeomvals = ionodegeom();
  }

  if(ionodegeomvals == NULL)return;

  /* Compute the number of nodes per I/O node along each direction */
  for(i = 0; i < 4; i++){
    if(dim_mach[i] % ionodegeomvals[i] != 0)status++;
    nodes_per_ionode[i] = dim_mach[i]/ionodegeomvals[i];
  }
  
  if(status){
    node0_printf("init_io_node: ionode geometry %d %d %d %d \n",
		 ionodegeomvals[0], ionodegeomvals[1],
		 ionodegeomvals[2], ionodegeomvals[3]);
    node0_printf("is incommensurate with node geometry %d %d %d %d\n",
		 dim_mach[0], dim_mach[1], dim_mach[2], dim_mach[3]);
    terminate(1);
  }
}

#endif

/*--------------------------------------------------------------------*/

/* Sets the QMP logical topology if we need one */
static void set_qmp_layout_grid(const int *geom, int n){
  if(geom == NULL)return;
  if(mynode()==0)printf("Setting QMP layout_grid to %d %d %d %d\n",
			geom[0], geom[1], geom[2], geom[3]);
  if(QMP_declare_logical_topology(geom, n) != QMP_SUCCESS){
    node0_printf("setup_layout: QMP_declare_logical_topology failed on %d %d %d %d \n",
		 geom[0], geom[1], geom[2], geom[3] );
    terminate(1);
  }
}

#if 0
/* Write my host name to a unique file for my node */
static void
mpi_whoami(void)
{
  char cmd[128];
  int m = mynode();
  
  sprintf(cmd,"hostname > node%d",m);
  system(cmd);

  sprintf(cmd,"echo ionode %d >> node%d",io_node(m),m);
  system(cmd);
}
#endif

/*------------------------------------------------------------------*/
/* Initialization entry point */

void
setup_layout(void)
{
  int c[4];
  int i,n_mach;
  int d[4];

#ifdef FIX_NODE_GEOM
  int *geom = node_geometry;
#else
  int *geom = NULL;
#endif

  if(mynode()==0){
    printf("LAYOUT = Hypercubes, options = ");
    printf("QDP");
    printf("\n");
  }

  /* Is there already a grid? 
     This could be a grid architecture with a preset dimension, or
     a geometry could have been set by the -qmp-geom command line arg. 
     In either case we have a nonzero allocated number of dimensions. 
  */

  if(QMP_get_allocated_number_of_dimensions() == 0)
    /* Set the geometry if requested */
    set_qmp_layout_grid(geom, 4);

  c[0] = nx;
  c[1] = ny;
  c[2] = nz;
  c[3] = nt;
  QDP_set_latsize(4, c);
  QDP_create_layout();
  sites_on_node = QDP_sites_on_node;
  even_sites_on_node = QDP_subset_len(QDP_even);
  odd_sites_on_node = QDP_subset_len(QDP_odd);
  n_mach = QMP_get_logical_number_of_dimensions();
  dim_mach = QMP_get_logical_dimensions();

  /* Initialize I/O node function */
#ifdef FIX_IONODE_GEOM
  init_io_node();
#endif
  
  /* Report sublattice dimensions */
  for(i = 0; i < 4; i++){
    /* Any extra machine dimensions are assumed to be 1 */
    if(i < n_mach)d[i] = c[i]/dim_mach[i];
    else d[i] = c[i];
  }
  if( mynode()==0)
    printf("ON EACH NODE %d x %d x %d x %d\n",d[0],d[1],d[2],d[3]);

#if 0
  mpi_whoami();  /* Debug */
#endif
}

int
node_number(int x, int y, int z, int t)
{
  int c[4];

  c[0] = x;
  c[1] = y;
  c[2] = z;
  c[3] = t;
  return QDP_node_number(c);
}

int
node_index(int x, int y, int z, int t)
{
  int c[4];

  c[0] = x;
  c[1] = y;
  c[2] = z;
  c[3] = t;
  return QDP_index(c);
}

size_t num_sites(int node) {
    return( sites_on_node );
}

const int *get_logical_dimensions(){
  return dim_mach;
}

const int *get_logical_coordinate(){
  return QMP_get_logical_coordinates();
}

/* Map node number and index to coordinates  */
void get_coords(int coords[], int node, int index){
  QDP_get_coords(coords, node, index);
}

/* io_node(node) maps a node to its I/O node.  The nodes are placed on
   a node lattice with dimensions dim_mach.  The I/O partitions are
   hypercubes of the node lattice.  The dimensions of the hypercube are
   given by nodes_per_ionode.  The I/O node is at the origin of that
   hypercube. */

/*------------------------------------------------------------------*/
/* Map any node to its I/O node */
int io_node(const int node){
  int i; 
  int io_node_coords[4];

  /* If we don't have I/O partitions, each node does its own I/O */
  if(ionodegeomvals == NULL)
    return node;

  /* Get the machine coordinates for the specified node */
  lex_coords(io_node_coords, 4, dim_mach, node);

  /* Round the node coordinates down to get the io_node coordinate */
  for(i = 0; i < 4; i++)
    io_node_coords[i] = nodes_per_ionode[i] * 
      (io_node_coords[i]/nodes_per_ionode[i]);
  
  /* Return the linearized machine coordinates of the I/O node */
  return (int)lex_rank(io_node_coords, 4, dim_mach);
}
