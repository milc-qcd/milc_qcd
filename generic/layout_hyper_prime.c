/******** layout_hyper_prime.c *********/
/* MIMD version 7 */
/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES */

/* The machine is viewed as a Cartesian grid of nodes.  Each node is,
   in turn, a Cartesian grid of MPI PEs (processor elements in MPI
   jargon), one rank per PE.  The combined subdivisions result
   in a overall Cartesian grid of MPI PEs.  Lattice sites are
   distributed uniformly in hypercubic sublattices across the MPI PEs.
   This file contains the procedures that define the mapping between a
   site coordinate and its MPI PE rank plus linear index within that
   rank.

   The division of the node sublattices and the further local
   subdivision within a node can be controlled by command-line
   parameters or by input parameters.  By default, the division is
   done automatically according to the hyper_prime algorithm by
   S. Gottlieb.

   NOTE: Some of the traditional terminology is apt to confuse.  The
   original versions of the code equated "node" with MPI rank.  For
   backward compatibility we have kept that terminology in global
   names.  See below.
   
   GLOBAL NAMES DEFINED HERE

   setup_layout() is the initialization call.  Determines which sites
                  go on which PE ranks and and sets the global
                  variables "sites_on_node", "even_sites_on_node",
                  "odd_sites_on_node"

   num_sites(rank) returns the number of sites on a PE rank

   node_number(x,y,z,t) returns the PE rank on which the site with
                        coordinates x,y,z,t lives.  Old-style name. 
			This procedure should be called "pe_rank_from_coords".

   node_index(x,y,z,t) returns the index of the site on its PE rank -
                       ie the site is lattice[node_index(x,y,z,t)].
                       Old-style name.  This procedure should be
                       called "site_index_from_coords"

   get_logical_dimensions() returns the dimensions of the PE grid.

   get_logical_coordinates() returns the PE grid coordinates of this
                           MPI rank

   get_coords() returns the site coordinates for a given PE rank and index
                         (the inverse of node_number + node_index)

   io_node(rank) maps a PE rank to its I/O rank (used if there are I/O
                    partitions).  Old-style name.  Should be called
		    "io_pe_rank".

   sites_on_node       Number of sites on this PE
   even_sites_on_node  Number of even sites on this PE
   odd_sites_on_node   Number of odd sites on this PE

*/

#include "generic_includes.h"
#ifdef HAVE_GRID
#include "../include/generic_grid.h"
#endif
#ifdef HAVE_QMP
#include <qmp.h>
#endif

static int squaresize[4];	   /* local dimensions (number of
				      sites) of the sublattice on one
				      node (hypercube) */
static int nsquares[4];	           /* number of site hypercubes in each direction, one hypercube per
				      node */
static int machine_coordinates[4]; /* logical machine (MPI PE rank) coordinates for this MPI rank */ 

static int pes_per_iopartition[4];    /* dimensions (number of MPI PEs) in one 
				      I/O partition. Must be a divisor of nsquares */
static int *ionodegeomvals = NULL; /* number of I/O partitions in each direction */

int prime[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};
# define MAXPRIMES ( sizeof(prime) / sizeof(int) )

/*------------------------------------------------------------------*/
/* Convert PE rank to coordinates */
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
/* Parity of the coordinate */
static int coord_parity(int r[]){
  return (r[0] + r[1] + r[2] + r[3]) % 2;
}

/*------------------------------------------------------------------*/
/* Convert coordinate to linear lexicographic PE rank (inverse of
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

#ifdef HAVE_QMP

/*--------------------------------------------------------------------*/
/* Sets the QMP logical topology from "geom" */
static void set_qmp_logical_topology(const int *geom, int n){

  /* Has a geometry already been specified by the -geom command-line
     argument or on the input parameter line "node_geometry"? */
  /* If not, don't set the grid geometry here */
  if(geom == NULL)return;

  /* If so, then pass the grid dimensions to QMP now */
  if(QMP_declare_logical_topology(geom, n) != QMP_SUCCESS){
    if(mynode()==0)printf("setup_layout: QMP_declare_logical_topology failed on %d %d %d %d \n",
		 geom[0], geom[1], geom[2], geom[3] );
    terminate(1);
  }
}

/*--------------------------------------------------------------------*/
static void setup_qmp_grid(const int *nsquares2, int ndim2){
  int ndim = 4;
  int len[4] = {nx, ny, nz, nt};
  int i;

  if(mynode()==0){
    printf("qmp_grid,");
    printf("\n");
  }

  for(i=0; i<ndim; i++) {
    if(i<ndim2) nsquares[i] = nsquares2[i];
    else nsquares[i] = 1;
  }

  if(mynode()==0){
    printf("Using machine geometry: ");
    for(i=0; i<ndim; i++){
      printf("%d ",nsquares[i]);
      if(i < ndim-1)printf("X ");
    }
    printf("\n");
  }

  for(i=0; i<ndim; i++) {
    if(len[i]%nsquares[i] != 0) {
      if(mynode()==0)printf("LATTICE SIZE DOESN'T FIT GRID\n");
      fflush(stdout);
      QMP_abort(0);
    }
    squaresize[i] = len[i]/nsquares[i];
  }
}
#endif

/*--------------------------------------------------------------------*/
/* This version divides the lattice by factors of prime numbers in any of the
   four directions.  It prefers to divide the longest dimensions,
   which mimimizes the area of the surfaces.  Similarly, it prefers
   to divide dimensions which have already been divided, thus not
   introducing more off-node directions.

	S. Gottlieb, May 18, 1999
	The code will start trying to divide with the largest prime factor
	and then work its way down to 2.  The current maximum prime is 53.
	The array of primes on line 46 may be extended if necessary.
*/

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
    
    /* find largest dimension of h-cubes divisible by prime[k] */
    for(j=0,dir=XUP;dir<=TUP;dir++)
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

void setup_fixed_geom(int const *geom, int n){
  int i;
  int node_count;
  int len[4];
  int status;
  char myname[] = "setup_fixed_geom";

// #ifdef FIX_NODE_GEOM
//   if(geom != NULL){
//     if(mynode()==0)printf("%s: Preallocated machine geometry overrides request\n", myname);
//   }
// #endif

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
    if(mynode()==0)printf("\n%s: Requested geometry %d %d %d %d ", myname,
		 geom[0], geom[1], geom[2], geom[3]);
    if(mynode()==0)printf("does not match number of nodes %d\n",numnodes());
    terminate(1);
  }

  if(status){
    if(mynode()==0)printf("%s: Requested geometry %d %d %d %d ", myname,
		 geom[0], geom[1], geom[2], geom[3]);
    if(mynode()==0)printf("is not commensurate with the lattice dims %d %d %d %d\n",
		 nx, ny, nz, nt);
    terminate(1);
  }
}

/*------------------------------------------------------------------*/
/* Initialize io_node function */

#ifdef FIX_IONODE_GEOM

static void init_io_node(){
  int i;
  int status = 0;
  char myname[] = "init_io_node";

  /* 1. If FIX_IONODE_GEOM is in force, use the ionode_geometry in the parameter
        input file.

     2. Otherwise, if -ionodes is specified on the command line, use it 

     3. Otherwise, all PEs write parallel I/O.
  */

  if(ionodegeom() == NULL){
    ionodegeomvals = ionode_geometry;
  } else {
    if(mynode()==0)printf("%s: Using command line -ionodes geometry\n", myname);
    ionodegeomvals = ionodegeom();
  }

  if(ionodegeomvals == NULL)return;

  /* Compute the dimensions (number of PE ranks) for one I/O
     partition along each direction */
  for(i = 0; i < 4; i++){
    if(nsquares[i] % ionodegeomvals[i] != 0)status++;
    pes_per_iopartition[i] = nsquares[i]/ionodegeomvals[i];
  }
  
  if(status){
    if(mynode()==0)printf("%s: ionode geometry %d %d %d %d \n", myname,
		 ionodegeomvals[0], ionodegeomvals[1],
		 ionodegeomvals[2], ionodegeomvals[3]);
    if(mynode()==0)printf("is incommensurate with node geometry %d %d %d %d\n",
		 nsquares[0], nsquares[1], nsquares[2], nsquares[3]);
    terminate(1);
  }
}
#endif

/*------------------------------------------------------------------*/
/* Sets the MPI PE layout nsquares, squaresize */
/* For QMP, declares logical "topology" */

static void set_topology(){
  int k = mynode();
  int nd = 0;
  int const *geom;

  if(k == 0) printf("LAYOUT = Hypercubes, options = ");

  /*--------------------------------*/
#ifdef HAVE_QMP

  /* QMP treatment */

  /* The PE layout dimensions (geometry) are set as follows:

     1. If -qmp-geom is specified on the command line, use the QMP allocated geometry

     2. Otherwise, if -qmp-alloc-map is specificed on the command
     line, use the allocated dimensions.

     3. Otherwise if FIX_NODE_GEOM is in force and node_geometry (from
        the parameter input) is defined, use node_geometry

     4. Otherwise use the layout_hyper_prime algorithm to set the geometry
  */
  
  nd = QMP_get_number_of_job_geometry_dimensions();
  if(nd > 0){
    /* Use QMP job geometry */
    geom = QMP_get_job_geometry();
    setup_qmp_grid(geom, nd);
    if(mynode()==0){
      printf("QMP using job_geometry_dimensions");
      for(int j = 0; j < nd; j++)
	printf(" %d",geom[j]);
      printf("\n");
      fflush(stdout);
    }
  } else {
    nd = QMP_get_allocated_number_of_dimensions();
    if(nd > 0) {
      geom = QMP_get_allocated_dimensions();
      /* use allocated geometry */
      setup_qmp_grid(geom, nd);
      if(mynode()==0)printf("QMP using allocated_dimension\n");
      fflush(stdout);
    } else {
#ifdef FIX_NODE_GEOM
      nd = 4;
      geom = node_geometry;
      /* take geometry from input parameter node_geometry line */
      setup_fixed_geom(geom, nd);
      if(mynode()==0)printf("QMP with input parameter node_geometry\n");
      fflush(stdout);
#else
      setup_hyper_prime();
      nd = 4;
      geom = nsquares;
      if(mynode()==0)printf("QMP with automatic hyper_prime layout\n");
      fflush(stdout);
#endif
    }
  }

  set_qmp_logical_topology(geom, nd);
  
  /*--------------------------------*/
#else
  
  /* Non QMP treatment */
  
  /* The layout dimensions (geometry) are set as follows:

     1. If the command line has -geom use that geometry

     2. Otherwise, if FIX_NODE_GEOM is in force and the node_geometry
     parameters are specified as input parameters, use them.

     3. Otherwise set the geometry with the layout_hyper_prime 
     algorithm
  */
  
  nd = 4;
  geom = nodegeom();  /* Command line values */
  
#ifdef FIX_NODE_GEOM
  if(geom == NULL){
    geom = node_geometry; /* Input parameter values */
  }
#endif
  
  if(geom != NULL){
    /* Set the sublattice dimensions according to the specified geometry */
    if(mynode()==0)printf("with fixed input-parameter node_geometry\n");
    setup_fixed_geom(geom, nd);
  } else {
    /* Set the sublattice dimensions according to the hyper_prime algorithm */
    setup_hyper_prime();
    if(mynode()==0)printf("automatic hyper_prime layout\n");
  }
  
#endif

}
  
/*------------------------------------------------------------------*/
/* Initialization entry point */

void setup_layout(){

  /* Set topology: nsquares, squaresize */

  set_topology();

#ifdef HAVE_GRID

  /* Initlalize Grid */
  initialize_grid();

  /* Grid assigns its own machine coordinates */

  /* Find my PE rank, according to Grid */
  setup_grid_communicator(nsquares);
  int *pePos = query_grid_node_mapping();
  int peRank = grid_rank_from_processor_coor(pePos[0], pePos[1], pePos[2], pePos[3]);

  /* Reassign my rank with the communicator */
  reset_machine_rank(peRank);

  set_topology();

#endif

  /* Initialize I/O node function */
#ifdef FIX_IONODE_GEOM
  init_io_node();
#endif
  
  /* Compute machine coordinates for this node */
#ifdef HAVE_GRID
  grid_coor_from_processor_rank(machine_coordinates, mynode());
#else
  lex_coords(machine_coordinates, 4, nsquares, mynode());
#endif

  /* Number of sites on node */
  sites_on_node =
    squaresize[XUP]*squaresize[YUP]*squaresize[ZUP]*squaresize[TUP];

  if( mynode()==0)
    printf("ON EACH NODE (RANK) %d x %d x %d x %d\n",squaresize[XUP],squaresize[YUP],
	   squaresize[ZUP],squaresize[TUP]);
  even_sites_on_node = sites_on_node/2;

  /* Allow an odd number of sites per node
     With an odd number, the number of even and 
     odd sites varies with each node...
     If the machine coordinates are even, we have an extra even site */
  if( (sites_on_node%2 != 0) && 
      ( coord_parity(machine_coordinates) == 0))
    even_sites_on_node++;

  odd_sites_on_node = sites_on_node - even_sites_on_node;

}

/*------------------------------------------------------------------*/
/* Return the MPI PE rank that contains the site x, y, z, t */

int node_number(int x, int y, int z, int t) {
  register int i;
  x /= squaresize[XUP]; y /= squaresize[YUP];
  z /= squaresize[ZUP]; t /= squaresize[TUP];
#ifdef HAVE_GRID
  i = grid_rank_from_processor_coor(x, y, z, t);
#else
#ifdef HAVE_QMP
  int proc_coords[4] = {x, y, z, t};
  i = QMP_get_node_number_from(proc_coords);
#else
  i = x + nsquares[XUP]*( y + nsquares[YUP]*( z + nsquares[ZUP]*( t )));
#endif
#endif
  return i;
}

/*------------------------------------------------------------------*/
/* Return the serialized site index for the site x, y, z, t on the PE
   that has it */

size_t node_index(int x, int y, int z, int t) {
    size_t i,xr,yr,zr,tr;
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
size_t get_even_sites_on_node(int node){
  int mc[4];
  size_t ir;
  size_t meo, neven;
  int k = node;
  
  /* mc = the machine coordinates for node k */
#ifdef HAVE_QMP
  QMP_comm_get_logical_coordinates_from2(QMP_comm_get_default(), mc, k);
#else
  lex_coords(mc, 4, nsquares, k);
#endif

  /* meo = the parity of the machine coordinate */
  meo = coord_parity(mc);

  /* neven = the number of even sites on node k */
  neven = (sites_on_node + 1 - meo)/2;

  return neven;
}  

/*------------------------------------------------------------------*/
/* Map PE rank number and serialize site index to coordinates */
/* (The inverse of node_number and node_index) */
/* Assumes even sites come first */
void get_coords(int coords[], int node, size_t index){
  int mc[4];
  size_t ir;
  size_t meo, neven, xeo;
  int k = node;

  /* mc = the machine coordinates for node k */
#ifdef HAVE_GRID
  grid_coor_from_processor_rank(mc, k);
#else
#ifdef HAVE_QMP
  QMP_comm_get_logical_coordinates_from2(QMP_comm_get_default(), mc, k);
#else
  lex_coords(mc, 4, nsquares, k);
#endif
#endif

  /* meo = the parity of the machine coordinate */
  meo = coord_parity(mc);

  /* neven = the number of even sites on node k */
  neven = (sites_on_node + 1 - meo)/2;
  
  /* ir = the even part of the lexicographic index within the
     sublattice on node k */
  if(index < neven){
    ir = 2*index;
    xeo = 0;
  } else {
    ir = 2*(index - neven);
    xeo = 1;
  }

  /* coords = the sublattice coordinate */
  lex_coords(coords, 4, squaresize, ir);

  /* Add offset to get full lattice coordinate (still a 2-fold ambiguity) */
  coords[XUP] += mc[XUP]*squaresize[XUP];
  coords[YUP] += mc[YUP]*squaresize[YUP];
  coords[ZUP] += mc[ZUP]*squaresize[ZUP];
  coords[TUP] += mc[TUP]*squaresize[TUP];

  /* Adjust coordinate according to parity */
  if( coord_parity(coords) != xeo ){
    coords[XUP]++;
    if(coords[XUP] >= squaresize[XUP]*(mc[XUP]+1)){
      coords[XUP] -= squaresize[XUP];
      coords[YUP]++;
      if(coords[YUP] >= squaresize[YUP]*(mc[YUP]+1)){
	coords[YUP] -= squaresize[YUP];
	coords[ZUP]++;
	if(coords[ZUP] >= squaresize[ZUP]*(mc[ZUP]+1)){
	  coords[ZUP] -= squaresize[ZUP];
	  coords[TUP]++;
	}
      }
    }
  }

  /* Consistency checks for debugging */
  if((k = node_number(coords[0], coords[1], coords[2], coords[3])) 
     != node){
    printf("get_coords: coords %d %d %d %d for node %d index %lu map to wrong node %d\n",
	   coords[0], coords[1], coords[2], coords[3], node, index, k);
    terminate(1);
  }
  size_t kk;
  if((kk = node_index(coords[0], coords[1], coords[2], coords[3]))
      != index){
    printf("get_coords: coords %d %d %d %d for node %d index %lu map to wrong index %lu\n",
	   coords[0], coords[1], coords[2], coords[3], node, index, kk);
    terminate(1);
  }
}

/*------------------------------------------------------------------*/
/* Map any PE rank to its I/O PE rank */

/* io_node(node) maps an MPI PE rank to its I/O PE rank.  The ranks
   are placed on a PE lattice with dimensions nsquares.  The I/O
   partitions are hypercubes of the PE lattice.  The dimensions of
   the hypercube are given by pes_per_iopartition.  The I/O PE is at
   the origin of that hypercube. */

int io_node(const int node){
  int i; 
  int io_node_coords[4];

  /* If we don't have I/O partitions, each MPI rank does its own I/O */
  if(ionodegeomvals == NULL)
    return node;

  /* Get the PE coordinates for the specified PE rank */
#ifdef HAVE_GRID
  grid_coor_from_processor_rank(io_node_coords, node);
#else
  lex_coords(io_node_coords, 4, nsquares, node);
#endif

  /* Round the PE coordinates down to get the io_node coordinate */
  for(i = 0; i < 4; i++)
    io_node_coords[i] = pes_per_iopartition[i] * 
      (io_node_coords[i]/pes_per_iopartition[i]);
  
  /* Return the linearized machine coordinates of the I/O node */
#ifdef HAVE_GRID
  return grid_rank_from_processor_coor(io_node_coords[0], 
	     io_node_coords[1], io_node_coords[2], io_node_coords[3]);
#else
  return (int)lex_rank(io_node_coords, 4, nsquares);
#endif
}


