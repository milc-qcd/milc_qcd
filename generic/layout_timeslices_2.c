/******** layout_timeslices_2.c *********/
/* MIMD version 7 */
/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES */
/* NOT IN THE CURRENT TEST STREAM.  USE WITH CAUTION! */

/* This version puts entire timeslices on nodes, if it can. It
   requires that nt, the time extent, is a multiple of the number of
   nodes used, or that the number of nodes be a multiple of nt, and
   that z be divisble by the factor..  We hope this speeds up spatial
   FFT's.

   Created May 28, 1998 S. Gottlieb

   3/29/00 EVENFIRST is the rule now. CD.
   11/16/06 Added get_logical_dimensions and get_logical_coordinates
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
   These routines will change as we change our minds about how to distribute
     sites among the nodes.  Hopefully the setup routines will work for any
     consistent choices. (ie node_index should return a different value for
     each site on the node.)
*/
#include "generic_includes.h"
#ifdef HAVE_QMP
#include <qio.h>
#endif

static int nsquares[4];	           /* number of hypercubes in each direction */
static int machine_coordinates[4]; /* logical machine coordinates */ 
static int ntslices;		/* number of timeslices per node */
static int nzslices;		/* number of z values per node */
static int zcuts;	/* number of times we must cut in z-direction, i.e,
							zcuts*nzslices=nz */

void setup_layout(){

#ifdef HAVE_QMP
 if(QMP_get_msg_passing_type()==QMP_GRID){
   printf("This layout should not be used on a grid architecture\n");
   terminate(1);
 }
#endif

    if(mynode()==0){
	printf("LAYOUT = Timeslices, options = ");
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
	ntslices = 1;
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
if( mynode()==0)
  printf("ON EACH NODE %d x %d x %d x %d\n",nx,ny,nzslices,ntslices);
if( mynode()==0 && sites_on_node%2 != 0)
	printf("WATCH OUT FOR EVEN/ODD SITES ON NODE BUG!!!\n");
   even_sites_on_node = odd_sites_on_node = sites_on_node/2;

    /* Define geometry in case someone asks */

    nsquares[XUP] = nsquares[YUP] = 1;
    nsquares[ZUP] = zcuts;
    nsquares[TUP] = numnodes()/zcuts;

    machine_coordinates[XUP] = 0;
    machine_coordinates[YUP] = 0;
    machine_coordinates[ZUP] = mynode() % nzcuts;
    machine_coordinates[TUP] = mynode()/zcuts;
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

size_t num_sites(int node) {
    return( sites_on_node );
}

const int *get_logical_dimensions(){
  return nsquares;
}

const int *get_logical_coordinate(){
  return machine_coordinates;
}
