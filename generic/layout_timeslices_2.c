/******** layout_timeslices_2.c *********/
/* MIMD version 6 */
/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES */

/* This version puts entire timeslices on nodes, if it can. It
   requires that nt, the time extent, is a multiple of the number of
   nodes used, or that the number of nodes be a multiple of nt, and
   that z be divisble by the factor..  We hope this speeds up spatial
   FFT's.

   Created May 28, 1998 S. Gottlieb

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

int ntslices;		/* number of timeslices per node */
int nzslices;		/* number of z values per node */
int zcuts;		/* number of times we must cut in z-direction, i.e,
							zcuts*nzslices=nz */

void setup_layout(){
    if(mynode()==0){
	printf("LAYOUT = Timeslices, options = ");
	printf("EVENFIRST,");
	printf("\n");
    }

    if( nt%numnodes() !=0 ){ /* then maybe numnodes is a multiple of nt */
      if (numnodes()%nt == 0 && nz%(numnodes()/nt) == 0 ) ;/* we can 
							      do the layout */
		elseif(mynode()==0)printf(
	    "LAYOUT: Can't lay out this lattice: nt not multiple of nummodes\n");
	    		terminate(1);
	}

    if( nt >= numnodes() ) {/* we only have to divide in t-direction */
    	ntslices = nt / numnodes();
	zcuts = 1;
    	nzslices = nz;
    }
    else {
	ntslizes = 1;
	zcuts = numnodes()/nt;
	nzslices = nz/zcuts;
    }
    sites_on_node = nx*ny*nzslices*ntslices;
    /* Need even number of sites per hypercube */
    if( mynode()==0)if( sites_on_node%2 != 0){
	printf("SORRY, CAN'T LAY OUT THIS LATTICE\n");
	terminate(0);
    }
    even_sites_on_node = odd_sites_on_node = sites_on_node/2;
}

int node_number(int x, int y, int z, int t) {
    t /= ntslices;
    z /= nzslices;
    return( t*zcuts + z );
}

int node_index(int x, int y, int z, int t) {
register int i, tr, zr;
    tr = t % ntslices;
    zr = z % nzslices;
    i = x + nx*(y + ny*(zr + nzslices*tr));
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


