/******** layout_qcdoc.c *********/
/* MIMD version 7 */
/* NO LONGER SUPPORTED C.D. 8/22/10*/
/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES */

/* This version divides each dimension by a machine determined
   grid size in that dimension.  Works only with QMP, because it
   uses a QMP call to get grid dimensions.
   Fails if lattice dimensions aren't divisible by grid dimensions.
   9/4/02 DT and EG
*/

/*
   setup_layout() does any initial setup.  When it is called the
     lattice dimensions nx,ny,nz and nt have been set.
     This routine sets the global variables "sites_on_node",
     "even_sites_on_node" and "odd_sites_on_node".
   num_sites(node) returns the number of sites on a node
   node_number(x,y,z,t) returns the node number on which a site lives.
   node_index(x,y,z,t) returns the index of the site on the node - ie the
     site is lattice[node_index(x,y,z,t)].
   get_logical_dimensions() returns the machine dimensions
   get_logical_coordinates() returns the mesh coordinates of this node
   get_coords() returns the coordinates for a given node and index
       (the inverse of node_number + node_index)
   io_node(node) maps nodes to their I/O node (for I/O partitions)
   These routines will change as we change our minds about how to distribute
     sites among the nodes.  Hopefully the setup routines will work for any
     consistent choices. (ie node_index should return a different value for
     each site on the node.)
*/
#include <qmp.h>
#include "generic_includes.h"

static int machine_dimensions[ 4 ];  /* DBR: Should this be 4? */
static int machine_nx = -1;
static int machine_ny = -1;
static int machine_nz = -1;
static int machine_nt = -1;

static int sub_lattice_dimensions[ 4 ];
static int sub_lattice_nx = -1;
static int sub_lattice_ny = -1;
static int sub_lattice_nz = -1;
static int sub_lattice_nt = -1;

static int sub_lattice_volume = -1;

/*------------------------------------------------------------------*/
/* Convert rank to coordinates */
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

#ifdef FIX_IONODE_GEOM

/*------------------------------------------------------------------*/
/* Convert coordinate to linear lexicographic rank (inverse of
   lex_coords) */

static size_t lex_rank(const int coords[], int dim, int size[])
{
  int d;
  size_t rank = coords[dim-1];

  for(d = dim-2; d >= 0; d--){
    rank = rank * size[d] + coords[d];
  }
  return rank;
}

#endif

#ifdef FIX_IONODE_GEOM

static int io_node_coords[4];
static int nodes_per_ionode[4];

/*------------------------------------------------------------------*/
/* Initialize io_node function */


static void init_io_node(){
  int i;
  int status = 0;

  /* Compute the number of nodes per I/O node along each direction */
  for(i = 0; i < 4; i++){
    if(machine_dimensions[i] % ionode_geometry[i] != 0)status++;
    nodes_per_ionode[i] = machine_dimensions[i]/ionode_geometry[i];
  }
  
  if(status){
    node0_printf("init_io_node: ionode geometry %d %d %d %d \n",
		 ionode_geometry[0], ionode_geometry[1],
		 ionode_geometry[2], ionode_geometry[3]);
    node0_printf("is incommensurate with node geometry %d %d %d %d\n",
		 machine_dimensions[0], machine_dimensions[1], machine_dimensions[3], machine_dimensions[3]);
    terminate(1);
  }
}
#endif

void setup_layout( void )
{
  const int *p_machine_dimensions = NULL;
  int number_machine_dimensions = -1;

  int i;

  number_machine_dimensions = QMP_get_logical_number_of_dimensions();
  printf( "number of QMP machine dimensions = %i\n", number_machine_dimensions );

  p_machine_dimensions = QMP_get_logical_dimensions();
  if( p_machine_dimensions == NULL )
  {
    printf( "p_machines_dimensions is NULL\n" );
    terminate( 0 );
  }

  for( i = 0; i < number_machine_dimensions; i++ )
  {
    printf( "QMP machine dimension ( %i ) = %i\n", i, p_machine_dimensions[ i ] );
  }

  machine_nx = p_machine_dimensions[ 0 ];
  machine_ny = p_machine_dimensions[ 1 ];
  machine_nz = p_machine_dimensions[ 2 ];
  machine_nt = p_machine_dimensions[ 3 ];

  p_machine_dimensions = NULL;

  machine_dimensions[ XUP ] = machine_nx;
  machine_dimensions[ YUP ] = machine_ny;
  machine_dimensions[ ZUP ] = machine_nz;
  machine_dimensions[ TUP ] = machine_nt;

  printf( "machine_nx = %i\n", machine_nx );
  printf( "machine_ny = %i\n", machine_ny );
  printf( "machine_nz = %i\n", machine_nz );
  printf( "machine_nt = %i\n", machine_nt );

  /* Each lattice dimension must be a mutliple of the corresponding machine dimension. */

  if( ( nx % machine_nx ) != 0 )
  {
    printf( "nx = %i is not a multiple of machine_nx = %i\n", nx, machine_nx );
    terminate( 0 );
  }
  if( ( ny % machine_ny ) != 0 )
  {
    printf( "ny = %i is not a multiple of machine_ny = %i\n", ny, machine_ny );
    terminate( 0 );
  }
  if( ( nz % machine_nz ) != 0 )
  {
    printf( "nz = %i is not a multiple of machine_nz = %i\n", nz, machine_nz );
    terminate( 0 );
  }
  if( ( nt % machine_nt ) != 0 )
  {
    printf( "nt = %i is not a multiple of machine_nt = %i\n", nt, machine_nt );
    terminate( 0 );
  }

  sub_lattice_nx = nx / machine_nx;
  sub_lattice_ny = ny / machine_ny;
  sub_lattice_nz = nz / machine_nz;
  sub_lattice_nt = nt / machine_nt;

  sub_lattice_dimensions[ XUP ] = sub_lattice_nx;
  sub_lattice_dimensions[ YUP ] = sub_lattice_ny;
  sub_lattice_dimensions[ ZUP ] = sub_lattice_nz;
  sub_lattice_dimensions[ TUP ] = sub_lattice_nt;

  printf( "sub_lattice_nx = %i\n", sub_lattice_nx );
  printf( "sub_lattice_ny = %i\n", sub_lattice_ny );
  printf( "sub_lattice_nz = %i\n", sub_lattice_nz );
  printf( "sub_lattice_nt = %i\n", sub_lattice_nt );

  sites_on_node = sub_lattice_nx * sub_lattice_ny * sub_lattice_nz * sub_lattice_nt;

  sub_lattice_volume = sites_on_node;

  /* The number of sites per node must be even. */

  if( mynode() == 0 )
  {
    if( sites_on_node % 2 != 0)
    {
        printf( "sites_on_node is not even\n" );
	terminate(0);
    }
  }

  even_sites_on_node = sites_on_node / 2;
  odd_sites_on_node  = sites_on_node / 2;

if( mynode()==0)
  printf("ON EACH NODE %d x %d x %d x %d\n",sub_lattice_nx,sub_lattice_ny,
                sub_lattice_nz,sub_lattice_nt);
if( mynode()==0 && sites_on_node%2 != 0)
	printf("WATCH OUT FOR EVEN/ODD SITES ON NODE BUG!!!\n");
}

int node_number( int x, int y, int z, int t )
{
  x /= sub_lattice_nx;  /* x -> machine_x */
  y /= sub_lattice_ny;  /* y -> machine_y */
  z /= sub_lattice_nz;  /* z -> machine_z */
  t /= sub_lattice_nt;  /* t -> machine_t */

  return( x + machine_nx * ( y + machine_ny * ( z + machine_nz * t ) ) );
}

int node_index( int x, int y, int z, int t )
{
  register int index;
  register int sub_lattice_x;
  register int sub_lattice_y;
  register int sub_lattice_z;
  register int sub_lattice_t;

  sub_lattice_x = x % sub_lattice_nx;
  sub_lattice_y = y % sub_lattice_ny;
  sub_lattice_z = z % sub_lattice_nz;
  sub_lattice_t = t % sub_lattice_nt;

  index = sub_lattice_x + sub_lattice_nx * ( sub_lattice_y + sub_lattice_ny * ( sub_lattice_z + sub_lattice_nz * sub_lattice_t ) );

  if( ( x + y + z + t ) % 2 == 0 )
  {
    /* even site */
    return( index / 2 );
  }
  else
  {
    /* odd site */
    return( ( index + sites_on_node ) / 2 );
  }
}

size_t num_sites( int node )
{
  return( sites_on_node );
}

const int *get_logical_dimensions(){
  return machine_dimensions;
}

const int *get_logical_coordinate(){
  return QMP_get_logical_coordinates();
}

/*------------------------------------------------------------------*/
/* Map node number and index to coordinates  */
void get_coords(int coords[], int node, int index){
  int mc[4];
  int ir;
  int eo;
  int k = node;

  /* Compute machine coordinates for node */
  lex_coords(mc, 4, machine_dimensions, k);

  /* Lexicographic index on node rounded to even */
  ir = 2*index;
  if(ir >= sites_on_node){
    ir -= sites_on_node;
    eo = 1;
  }
  else
    eo = 0;

  /* Convert to coordinates - result is two-fold ambiguous */
  lex_coords(coords, 4, sub_lattice_dimensions, ir);

  /* Adjust coordinate according to parity (assumes even sites_on_node) */
  if( (coords[0] + coords[1] + coords[2] + coords[3]) % 2 != eo){
    coords[0]++;
    if(coords[0] >= sub_lattice_dimensions[0]){
      coords[0] -= sub_lattice_dimensions[0];
      coords[1]++;
      if(coords[1] >= sub_lattice_dimensions[1]){
	coords[1] -= sub_lattice_dimensions[1];
	coords[2]++;
	if(coords[2] >= sub_lattice_dimensions[2]){
	  coords[2] -= sub_lattice_dimensions[2];
	  coords[3]++;
	}
      }
    }
  }

  /* Add offset for hypercube */
  coords[0] += mc[0]*sub_lattice_dimensions[0];
  coords[1] += mc[1]*sub_lattice_dimensions[1];
  coords[2] += mc[2]*sub_lattice_dimensions[2];
  coords[3] += mc[3]*sub_lattice_dimensions[3];

  /* Consistency checks for debugging */
  if((k = node_number(coords[0], coords[1], coords[2], coords[3])) 
     != node){
    printf("get_coords: coords %d %d %d %d for node %d map to wrong node %d\n",
	   coords[0], coords[1], coords[2], coords[3], node, k);
    terminate(1);
  }
  if((k = node_index(coords[0], coords[1], coords[2], coords[3]))
      != index){
    printf("get_coords: coords %d %d %d %d for index %d map to wrong index %d\n",
	   coords[0], coords[1], coords[2], coords[3], index, k);
    terminate(1);
  }
}

/* io_node(node) maps a node to its I/O node.  The nodes are placed on
   a node lattice with dimensions machine_dimensions.  The I/O partitions are
   hypercubes of the node lattice.  The dimensions of the hypercube are
   given by nodes_per_ionode.  The I/O node is at the origin of that
   hypercube. */

#ifdef FIX_IONODE_GEOM

/*------------------------------------------------------------------*/
/* Map any node to its I/O node */
int io_node(const int node){
  int i,j,k; 

  /* Get the machine coordinates for the specified node */
  lex_coords(io_node_coords, 4, machine_dimensions, node);

  /* Round the node coordinates down to get the io_node coordinate */
  for(i = 0; i < 4; i++)
    io_node_coords[i] = nodes_per_ionode[i] * 
      (io_node_coords[i]/nodes_per_ionode[i]);
  
  /* Return the linearized machine coordinates of the I/O node */
  return (int)lex_rank(io_node_coords, 4, machine_dimensions);
}

#else

/*------------------------------------------------------------------*/
/* If we don't have I/O partitions, each node does its own I/O */
int io_node(int node){
  return node;
}
#endif
