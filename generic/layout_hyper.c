/******** layout_hyper.c *********/
/* MIMD version 6 */
/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES */

/* This version divides the lattice by factors of two in any of the
   four directions.  It prefers to divide the longest dimensions,
   which mimimizes the area of the surfaces.  Similarly, it prefers
   to divide dimensions which have already been divided, thus not
   introducing more off-node directions.

   This requires that the lattice volume be divisible by the number
   of nodes, which is a power of two.

   3/29/00 EVENFIRST is the rule now. CD.

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

int squaresize[4];	/* dimensions of hypercubes */
int nsquares[4];	/* number of hypercubes in each direction */

void setup_layout(){
register int i,j,dir;
    if(mynode()==0){
	printf("LAYOUT = Hypercubes, options = ");
	printf("EVENFIRST,");
	printf("\n");
    }

    /* Figure out dimensions of rectangle */
    squaresize[XUP] = nx; squaresize[YUP] = ny;
    squaresize[ZUP] = nz; squaresize[TUP] = nt;
    nsquares[XUP] = nsquares[YUP] = nsquares[ZUP] = nsquares[TUP] = 1;

    i = 1;	/* current number of hypercubes */
    while(i<numnodes()){
	/* figure out which direction to divide */

	/* find largest even dimension of h-cubes */
	for(j=1,dir=XUP;dir<=TUP;dir++)
	    if( squaresize[dir]>j && squaresize[dir]%2==0 ) j=squaresize[dir];

	/* if one direction with largest dimension has already been
	   divided, divide it again.  Otherwise divide first direction
	   with largest dimension. */
	for(dir=XUP;dir<=TUP;dir++)
	    if( squaresize[dir]==j && nsquares[dir]>1 )break;
	if( dir > TUP)for(dir=XUP;dir<=TUP;dir++)
	    if( squaresize[dir]==j )break;
	/* This can fail if I run out of factors of 2 in the dimensions */
	if(dir > TUP){
	    if(mynode()==0)printf(
	    "LAYOUT: Can't lay out this lattice, not enough factors of 2\n");
	    terminate(1);
	}

	/* do the surgery */
	i*=2; squaresize[dir] /= 2; nsquares[dir] *= 2;
    }

    sites_on_node =
	    squaresize[XUP]*squaresize[YUP]*squaresize[ZUP]*squaresize[TUP];
    /* Need even number of sites per hypercube */
    if( mynode()==0)if( sites_on_node%2 != 0){
	printf("SORRY, CAN'T LAY OUT THIS LATTICE\n");
	terminate(0);
    }
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

int num_sites(int node) {
    return( sites_on_node );
}
