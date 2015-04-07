#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

#include "params.h"
#include "../include/random.h"   /* For double_prn */

typedef struct {
    /* The first part is standard to all programs */
	/* coordinates of this site */
	short x,y,z,t;
	/* is it even or odd? */
	char parity;
	/* my index in the array */
	int index;
	/* The state information for a random number generator */
  	double_prn site_prn;
	/* align to double word boundary (kludge for Intel compiler) */
	int space1;
} site;

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  int iseed;              /* Not used in this application */
EXTERN  int volume;		/* volume of lattice = nx*ny*nz*nt */

#ifdef FIX_NODE_GEOM
EXTERN  int node_geometry[4];  /* Specifies fixed "nsquares" (i.e. 4D
			    hypercubes) for the compute nodes in each
			    coordinate direction.  Must be divisors of
			    the lattice dimensions */
#ifdef FIX_IONODE_GEOM
EXTERN int ionode_geometry[4]; /* Specifies fixed "nsquares" for I/O
			     partitions in each coordinate direction,
			     one I/O node for each square.  The I/O
			     node is at the origin of the square.
			     Must be divisors of the node_geometry. */
#endif
#endif
EXTERN  params param;           /* user input parameters */

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

EXTERN char utc_date_time[64];
EXTERN char hostname[128];

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 16	/* Number of generic pointers */
/* NEED 8 WHEN GAUGEFIXING */
EXTERN char ** gen_pt[N_POINTERS];

EXTERN site *lattice;

#endif /* _LATTICE_H */
