/******** layout_hyper_tstretch.c *********/
/* MIMD version 7 */
/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES */
/* NOT IN THE CURRENT TEST STREAM.  USE WITH CAUTION! */

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

   DT 1/7/03 Allow last subvolume in time direction to have one fewer
   time slice than others.  For example, put a 24^3x64 lattice on 120
   processors by dividing time direction into five parts as 64 =
   4*13+12.  This is NOT implemented in the most general way.

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
   These routines will change as we change our minds about how to distribute
     sites among the nodes.  Hopefully the setup routines will work for any
     consistent choices. (ie node_index should return a different value for
     each site on the node.)
*/
#include "generic_includes.h"
#ifdef HAVE_QMP
#include <qio.h>
#endif

int squaresize[4];	/* dimensions of hypercubes */
int nsquares[4];	/* number of hypercubes in each direction */
static int sites_on_lastnodes;

int prime[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};
# define MAXPRIMES ( sizeof(prime) / sizeof(int) )

void setup_layout(){
register int i,j,k,dir;

#ifdef HAVE_QMP
 if(QMP_get_msg_passing_type()==QMP_GRID){
   printf("This layout should not be used on a grid architecture\n");
   terminate(1);
 }
#endif

    if(mynode()==0){
	printf("LAYOUT = Hypercubes, options = ");
	printf("tstretch,");
	printf("\n");
    }

    /* Figure out dimensions of rectangle */
    if( (nx*ny*nz*nt)%numnodes() == 0){
      squaresize[XUP] = nx; squaresize[YUP] = ny;
      squaresize[ZUP] = nz; squaresize[TUP] = nt;
    }
    else if ( (nx*ny*nz*(nt+1))%numnodes() == 0){
      squaresize[XUP] = nx; squaresize[YUP] = ny;
      squaresize[ZUP] = nz; squaresize[TUP] = nt+1;
    }
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
	    if(mynode()==0)printf(
	    "LAYOUT: Can't lay out this lattice, not enough factors of %d\n"
		,prime[k]);
	    terminate(1);
	}

	/* do the surgery */
	i*=prime[k]; squaresize[dir] /= prime[k]; nsquares[dir] *= prime[k];
    }

    sites_on_node =
	    squaresize[XUP]*squaresize[YUP]*squaresize[ZUP]*squaresize[TUP];
    sites_on_lastnodes = sites_on_node;

if( mynode()==0)
  printf("ON MOST NODES %d x %d x %d x %d\n",squaresize[XUP],squaresize[YUP],
                squaresize[ZUP],squaresize[TUP]);
    if( !((nx*ny*nz*nt)%numnodes() == 0) && (nx*ny*nz*(nt+1))%numnodes() == 0 ){
	/* stretched t direction by one */
	if( mynode()==0)printf("SOME NODES HAVE FEWER SITES\n");
	if( mynode() >=  nsquares[XUP]*nsquares[YUP]*nsquares[ZUP]*(nsquares[TUP]-1) ){
	    sites_on_node = squaresize[XUP]*squaresize[YUP]*squaresize[ZUP]*(squaresize[TUP]-1);
	    sites_on_lastnodes = sites_on_node;
	}
    }

    /* Need even number of sites per hypercube */
    if( mynode()==0)if( sites_on_node%2 != 0){
	printf("SORRY, CAN'T LAY OUT THIS LATTICE\n");
	terminate(0);
    }
if( mynode()==0 && sites_on_node%2 != 0)
	printf("WATCH OUT FOR EVEN/ODD SITES ON NODE BUG!!!\n");
    even_sites_on_node = odd_sites_on_node = sites_on_node/2;
}

int node_number(int x,int y,int z,int t) {
register int i;
    x /= squaresize[XUP]; y /= squaresize[YUP];
    z /= squaresize[ZUP]; t /= squaresize[TUP];
    i = x + nsquares[XUP]*( y + nsquares[YUP]*( z + nsquares[ZUP]*( t )));
    return( i );
}

int node_index(int x,int y,int z,int t) {
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

size_t num_sites(int node) {
  if( node >=  nsquares[XUP]*nsquares[YUP]*nsquares[ZUP]*(nsquares[TUP]-1) )
    return sites_on_lastnodes;
  else
    return sites_on_node;
}

