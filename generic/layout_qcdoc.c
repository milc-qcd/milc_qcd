/******** layout_qcdocsim.c *********/
/* MIMD version 6 */
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
   These routines will change as we change our minds about how to distribute
     sites among the nodes.  Hopefully the setup routines will work for any
     consistent choices. (ie node_index should return a different value for
     each site on the node.)
*/
#include "generic_includes.h"
#include <qmp.h>

int squaresize[4];	/* dimensions of hypercubes */
int nsquares[4];	/* number of hypercubes in each direction */
int *nsquares2;   /* HACK. Use QMP to get machine size */
int latdims[4];		/* size of lattice */
/*int *nsquares;*/

void setup_layout(){
register int i,j,dir;
#if 0
QMP_u32_t *QMP_get_allocated_dimensions();
#endif


/* nsquares=(int*)malloc(sizeof(int)*4);*/

#ifndef PPC440QCDOC
    if(mynode()==0){
	printf("LAYOUT = Grid, options = ");
	printf("EVENFIRST,");
	printf("\n");
    }
#endif /* PPC440QCDOC */

    printf("QMP num dim= %d, this_node=%d\n", QMP_get_allocated_number_of_dimensions (), this_node);fflush(stdout);

    /* Figure out dimensions of rectangle */
#ifndef DUMMYTEST
    /*
        if( QMP_get_allocated_number_of_dimensions () != 4){
      node0_printf("BONEHEAD! THIS WORKS ONLY ON 4D GRID MACHINE\n"); terminate(0);
        }
    */
    nsquares2 = QMP_get_allocated_dimensions ();
    nsquares[XUP]=nsquares2[0];
    nsquares[YUP]=nsquares2[1];
    nsquares[ZUP]=nsquares2[2];
    nsquares[TUP]=nsquares2[3];
#ifndef PPC440QCDOC
printf("nsquares[XUP]=%d\n",nsquares[XUP]);
printf("nsquares[YUP]=%d\n",nsquares[YUP]);
printf("nsquares[ZUP]=%d\n",nsquares[ZUP]);
printf("nsquares[TUP]=%d\n",nsquares[TUP]);
#endif /* PPC440QCDOC */
#else  /*DUMMYTEST*/
    nsquares[XUP]=2;
    nsquares[YUP]=1;
    nsquares[ZUP]=3;
    nsquares[TUP]=numnodes()/( 2*1*3 );
#endif /*DUMMYTEST*/

    /* set grid dimensions */
    /* find dimensions of each block */
    latdims[XUP]=nx; latdims[YUP]=ny; latdims[ZUP]=nz; latdims[TUP]=nt;
    for(dir=XUP;dir<=TUP;dir++){
	if( latdims[dir]%nsquares[dir] != 0){
	    node0_printf("LATTICE SIZE DOESN'T FIT GRID\n"); terminate(0);
	}
	squaresize[dir] = latdims[dir]/nsquares[dir];
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
