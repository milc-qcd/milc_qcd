/******** layout_hyper.c *********/
/* MIMD code version 3 */
/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES */
/* Version from Craig McNeile, obtained Apr 1996 */

/* This version divides the lattice by factors of two in any of the
   four directions.  It prefers to divide the longest dimensions,
   which mimimizes the area of the surfaces.  Similarly, it prefers
   to divide dimensions which have already been divided, thus not
   introducing more off-node directions.

   This requires that the lattice volume be divisible by the number
   of nodes, which is a power of two.

   With the "ACCORDION" option, the numbers of the cubes in the lattice
   are Gray coded.  THIS OPTION IS IGNORED FOR NOW.

   With the "GRAYCODE" option the node numbers are gray coded so that
   adjacent lattice regions will physically be on adjacent nodes
   in a hypercube architecture

   With the "EVENFIRST" option the even sites are listed contiguously
   in the first part of lattice[], and the odd sites in the last part.
   In this version there must be an equal number of even and odd sites
   in each cube - in other words one of the dimensions of the h-cube must
   be even.
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
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <su3.h>
#include LATDEF
#include <comdefs.h>

int squaresize[4];	/* dimensions of hypercubes */
int nsquares[4];	/* number of hypercubes in each direction */

void setup_layout(){
register int i,j,k,dir;
    if(mynode()==0){
	printf("LAYOUT = Hypercubes, options = ");
#ifdef ACCORDION
	printf("ACCORDION(ignored),");
#endif
#ifdef GRAYCODE
	printf("GRAYCODE,");
#endif
#ifdef EVENFIRST
	printf("EVENFIRST,");
#endif
	printf("\n");
    }

/******set up the layout_flag for quark propagator IO routines *****/

#ifdef GRAYCODE
  #ifdef EVENFIRST
    layout_flag = HYPER_GRAY_EVENFIRST ;
/* BJ 2/3/97 
   #else 
      #error layout_hyper code needs to be updated 
    #endif
  #else
    #error layout_hyper code needs to be updated  
*/ 
#endif
#endif

/****  ------------------------------  **********/


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
	    exit(1);
	}

	/* do the surgery */
	i*=2; squaresize[dir] /= 2; nsquares[dir] *= 2;
    }

    sites_on_node =
	    squaresize[XUP]*squaresize[YUP]*squaresize[ZUP]*squaresize[TUP];
#ifdef EVENFIRST
    /* Need even number of sites per hypercube */
    if( mynode()==0)if( sites_on_node%2 != 0){
	printf("SORRY, CAN'T LAY OUT THIS LATTICE\n");
	terminate(0);
    }
#endif
    even_sites_on_node = odd_sites_on_node = sites_on_node/2;
}

int node_number(x,y,z,t) int x,y,z,t; {
register int i;
    x /= squaresize[XUP]; y /= squaresize[YUP];
    z /= squaresize[ZUP]; t /= squaresize[TUP];
    i = x + nsquares[XUP]*( y + nsquares[YUP]*( z + nsquares[ZUP]*( t )));
#ifdef ACCORDION
    /* Not yet implemented */
#endif
#ifdef GRAYCODE
    return( i ^ (i>>1) );	/* Gray code of i */
#else
    return( i );
#endif
}

int node_index(x,y,z,t) int x,y,z,t; {
register int i,xr,yr,zr,tr;
    xr = x%squaresize[XUP]; yr = y%squaresize[YUP];
    zr = z%squaresize[ZUP]; tr = t%squaresize[TUP];
    i = xr + squaresize[XUP]*( yr + squaresize[YUP]*( zr + squaresize[ZUP]*tr));
#ifdef EVENFIRST
    if( (x+y+z+t)%2==0 ){	/* even site */
	return( i/2 );
    }
    else {
	return( (i + sites_on_node)/2 );
    }
#else
    return( i );
#endif
}

int num_sites(node) int node; {
    return( sites_on_node );
}

