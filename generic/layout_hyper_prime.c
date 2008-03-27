/******** layout_hyper_prime.c *********/
/* MIMD version 7 */
/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES */

/* This version divides the lattice by factors of prime numbers in any of the
   four directions.  It prefers to divide the longest dimensions,
   which mimimizes the area of the surfaces.  Similarly, it prefers
   to divide dimensions which have already been divided, thus not
   introducing more off-node directions.

	S. Gottlieb, May 18, 1999
	The code will start trying to divide with the largest prime factor
	and then work its way down to 2.  The current maximum prime is 53.
	The array of primes on line 46 may be extended if necessary.

   This requires that the lattice volume be divisible by the number
   of nodes.  Each dimension must be divisible by a suitable factor
   such that the product of the four factors is the number of nodes.

   3/29/00 EVENFIRST is the rule now. CD.
   12/10/00 Fixed so k = MAXPRIMES-1 DT
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
#include "generic_includes.h"
#ifdef HAVE_QMP
#include <qmp.h>
#endif

static int squaresize[4];	   /* dimensions of hypercubes */
static int nsquares[4];	           /* number of hypercubes in each direction */
static int machine_coordinates[4]; /* logical machine coordinates */ 

int prime[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};
# define MAXPRIMES ( sizeof(prime) / sizeof(int) )

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

#ifdef HAVE_QMP

/*--------------------------------------------------------------------*/
/* Sets the QMP logical topology if we need one */
static void set_qmp_layout_grid(const int *geom, int n){
  if(geom == NULL)return;
  if(QMP_declare_logical_topology(geom, n) != QMP_SUCCESS){
    node0_printf("setup_layout: QMP_declare_logical_topology failed on %d %d %d %d \n",
		 geom[0], geom[1], geom[2], geom[3] );
    terminate(1);
  }
}

/*--------------------------------------------------------------------*/
static void setup_qmp_grid(){
  int ndim = 4;
  int len[4];
  int ndim2, i;
  const int *nsquares2;

  len[0] = nx; len[1] = ny; len[2] = nz; len[3] = nt;

  if(mynode()==0){
    printf("qmp_grid,");
    printf("\n");
  }

  ndim2 = QMP_get_allocated_number_of_dimensions();
  nsquares2 = QMP_get_allocated_dimensions();

  /* If the dimensions are not already allocated, use the
     node_geometry request.  Otherwise a hardware or command line
     specification trumps the parameter input. */
#ifdef FIX_NODE_GEOM
  if(ndim2 == 0){
    ndim2 = 4;
    nsquares2 = node_geometry;
  }
  else{
    node0_printf("setup_qmp_grid: Preallocated machine geometry overrides request\n");
  }
#endif

  if(mynode()==0){
    printf("Using machine geometry: ");
    for(i=0; i<ndim; i++){
      printf("%d ",nsquares2[i]);
      if(i < ndim-1)printf("X ");
    }
    printf("\n");
  }

  /* In principle, we could now rotate coordinate axes */
  /* Save this for a future upgrade */

  set_qmp_layout_grid(nsquares2, ndim2);

  ndim2 = QMP_get_logical_number_of_dimensions();
  nsquares2 = QMP_get_logical_dimensions();

  for(i=0; i<ndim; i++) {
    if(i<ndim2) nsquares[i] = nsquares2[i];
    else nsquares[i] = 1;
  }

  for(i=0; i<ndim; i++) {
    if(len[i]%nsquares[i] != 0) {
      node0_printf("LATTICE SIZE DOESN'T FIT GRID\n");
      QMP_abort(0);
    }
    squaresize[i] = len[i]/nsquares[i];
  }
}
#endif

/*--------------------------------------------------------------------*/
static void setup_hyper_prime(){
  int i,j,k,dir;

  if(mynode()==0){
    printf("hyper_prime,");
    printf("\n");
  }

  /* Figure out dimensions of rectangle */
  squaresize[XUP] = nx; squaresize[YUP] = ny;
  squaresize[ZUP] = nz; squaresize[TUP] = nt;
  nsquares[XUP] = nsquares[YUP] = nsquares[ZUP] = nsquares[TUP] = 1;
  
  i = 1;	/* current number of hypercubes */
  while(i<numnodes()){
    /* figure out which prime to divide by starting with largest */
    k = MAXPRIMES-1;
    while( (numnodes()/i)%prime[k] != 0 && k>0 ) --k;
    /* figure out which direction to divide */
    
    /* find largest even dimension of h-cubes */
    for(j=1,dir=XUP;dir<=TUP;dir++)
      if( squaresize[dir]>j && squaresize[dir]%prime[k]==0 )
	j=squaresize[dir];
    
    /* if one direction with largest dimension has already been
       divided, divide it again.  Otherwise divide first direction
       with largest dimension. */
    for(dir=XUP;dir<=TUP;dir++)
      if( squaresize[dir]==j && nsquares[dir]>1 )break;
    if( dir > TUP)for(dir=XUP;dir<=TUP;dir++)
      if( squaresize[dir]==j )break;
    /* This can fail if I run out of prime factors in the dimensions */
    if(dir > TUP){
      if(mynode()==0)
	printf("LAYOUT: Can't lay out this lattice, not enough factors of %d\n"
	       ,prime[k]);
      terminate(1);
    }
    
    /* do the surgery */
    i*=prime[k]; squaresize[dir] /= prime[k]; nsquares[dir] *= prime[k];
  }
}

/*--------------------------------------------------------------------*/

void setup_fixed_geom(int *geom, int n){
  int i;
  int node_count;
  int len[4];
  int status;

  len[0] = nx; len[1] = ny; len[2] = nz; len[3] = nt;

  node_count = 1;
  status = 0;
  for(i = 0; i < 4; i++){
    nsquares[i] = geom[i];
    node_count *= geom[i];
    if(len[i] % nsquares[i] != 0)status++;
    squaresize[i] = len[i]/nsquares[i];
  }

  if(node_count != numnodes()){
    node0_printf("/nsetup_fixed_geom: Requested geometry %d %d %d %d ",
		 geom[0], geom[1], geom[2], geom[3]);
    node0_printf("does not match number of nodes %d\n",numnodes());
    terminate(1);
  }

  if(status){
    node0_printf("setup_fixed_geom: Requested geometry %d %d %d %d ",
		 geom[0], geom[1], geom[2], geom[3]);
    node0_printf("is not commensurate with the lattice dims %d %d %d %d\n",
		 nx, ny, nz, nt);
    terminate(1);
  }
}

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
    if(nsquares[i] % ionode_geometry[i] != 0)status++;
    nodes_per_ionode[i] = nsquares[i]/ionode_geometry[i];
  }
  
  if(status){
    node0_printf("init_io_node: ionode geometry %d %d %d %d \n",
		 ionode_geometry[0], ionode_geometry[1],
		 ionode_geometry[2], ionode_geometry[3]);
    node0_printf("is incommensurate with node geometry %d %d %d %d\n",
		 nsquares[0], nsquares[1], nsquares[3], nsquares[3]);
    terminate(1);
  }
}
#endif

/*------------------------------------------------------------------*/
/* Initialization entry point */

void setup_layout(){
  int k = mynode();
#ifdef FIX_NODE_GEOM
  int *geom = node_geometry;
#else
  int *geom = NULL;
#endif

  if(k == 0)
    printf("LAYOUT = Hypercubes, options = ");

#ifdef HAVE_QMP
  /* QMP treatment */
  /* Is there already a grid? 
     This could be a grid architecture with a preset dimension, or
     a geometry could have been set by the -qmp-geom command line arg. 
     In either case we have a nonzero allocated number of dimensions. 
*/
  if(QMP_get_allocated_number_of_dimensions() == 0)
    /* Set the geometry if requested */
    set_qmp_layout_grid(geom, 4);

  /* Has a grid been set up now? */
  if(QMP_get_msg_passing_type() == QMP_GRID)
    setup_qmp_grid();
  else if(geom != NULL)
    setup_fixed_geom(geom, 4);
  else
    setup_hyper_prime();

#else

  /* Non QMP treatment */
  if(geom != NULL)
    setup_fixed_geom(geom, 4);
  else
    setup_hyper_prime();

#endif

#ifdef FIX_IONODE_GEOM
  /* Initialize I/O node function */
  init_io_node();
#endif
  
  /* Compute machine coordinates for this node */
  lex_coords(machine_coordinates, 4, nsquares, k);

  /* Number of sites on node */
  sites_on_node =
    squaresize[XUP]*squaresize[YUP]*squaresize[ZUP]*squaresize[TUP];
  /* Need even number of sites per hypercube */
  if( mynode()==0)if( sites_on_node%2 != 0){
    printf("SORRY, CAN'T LAY OUT THIS LATTICE\n");
    terminate(0);
  }
  if( mynode()==0)
    printf("ON EACH NODE %d x %d x %d x %d\n",squaresize[XUP],squaresize[YUP],
	   squaresize[ZUP],squaresize[TUP]);
  if( mynode()==0 && sites_on_node%2 != 0)
    printf("WATCH OUT FOR EVEN/ODD SITES ON NODE BUG!!!\n");
  even_sites_on_node = odd_sites_on_node = sites_on_node/2;
}

/*------------------------------------------------------------------*/
int node_number(int x, int y, int z, int t) {
register int i;
    x /= squaresize[XUP]; y /= squaresize[YUP];
    z /= squaresize[ZUP]; t /= squaresize[TUP];
    i = x + nsquares[XUP]*( y + nsquares[YUP]*( z + nsquares[ZUP]*( t )));
    return( i );
}

/*------------------------------------------------------------------*/
int node_index(int x, int y, int z, int t) {
register int i,xr,yr,zr,tr;
    xr = x%squaresize[XUP]; yr = y%squaresize[YUP];
    zr = z%squaresize[ZUP]; tr = t%squaresize[TUP];
    i = xr + squaresize[XUP]*( yr + squaresize[YUP]*( zr + squaresize[ZUP]*tr));
    if( (x+y+z+t)%2==0 ){	/* even site */
	return( i/2 );
    }
    else {
	return( (i + sites_on_node)/2 );
    }
}

/*------------------------------------------------------------------*/
size_t num_sites(int node) {
    return( sites_on_node );
}

/*------------------------------------------------------------------*/
const int *get_logical_dimensions(){
  return nsquares;
}

/*------------------------------------------------------------------*/
/* Coordinates simulate a mesh architecture and must correspond
   to the node_number result */
const int *get_logical_coordinate(){
  return machine_coordinates;
}

/*------------------------------------------------------------------*/
/* Map node number and index to coordinates  */
void get_coords(int coords[], int node, int index){
  int mc[4];
  int ir;
  int eo;
  int k = node;

  /* Compute machine coordinates for node */
  lex_coords(mc, 4, nsquares, k);

  /* Lexicographic index on node rounded to even */
  ir = 2*index;
  if(ir >= sites_on_node){
    ir -= sites_on_node;
    eo = 1;
  }
  else
    eo = 0;

  /* Convert to coordinates - result is two-fold ambiguous */
  lex_coords(coords, 4, squaresize, ir);

  /* Adjust coordinate according to parity (assumes even sites_on_node) */
  if( (coords[XUP] + coords[YUP] + coords[ZUP] + coords[TUP]) % 2 != eo){
    coords[XUP]++;
    if(coords[XUP] >= squaresize[XUP]){
      coords[XUP] -= squaresize[XUP];
      coords[YUP]++;
      if(coords[YUP] >= squaresize[YUP]){
	coords[YUP] -= squaresize[YUP];
	coords[ZUP]++;
	if(coords[ZUP] >= squaresize[ZUP]){
	  coords[ZUP] -= squaresize[ZUP];
	  coords[TUP]++;
	}
      }
    }
  }

  /* Add offset for hypercube */
  coords[XUP] += mc[XUP]*squaresize[XUP];
  coords[YUP] += mc[YUP]*squaresize[YUP];
  coords[ZUP] += mc[ZUP]*squaresize[ZUP];
  coords[TUP] += mc[TUP]*squaresize[TUP];

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
   a node lattice with dimensions nsquares.  The I/O partitions are
   hypercubes of the node lattice.  The dimensions of the hypercube are
   given by nodes_per_ionode.  The I/O node is at the origin of that
   hypercube. */

#ifdef FIX_IONODE_GEOM

/*------------------------------------------------------------------*/
/* Map any node to its I/O node */
int io_node(const int node){
  int i,j,k; 

  /* Get the machine coordinates for the specified node */
  lex_coords(io_node_coords, 4, nsquares, node);

  /* Round the node coordinates down to get the io_node coordinate */
  for(i = 0; i < 4; i++)
    io_node_coords[i] = nodes_per_ionode[i] * 
      (io_node_coords[i]/nodes_per_ionode[i]);
  
  /* Return the linearized machine coordinates of the I/O node */
  return (int)lex_rank(io_node_coords, 4, nsquares);
}

#else

/*------------------------------------------------------------------*/
/* If we don't have I/O partitions, each node does its own I/O */
int io_node(int node){
  return node;
}
#endif
