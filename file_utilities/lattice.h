#ifndef _LATTICE_H
#define _LATTICE_H
/* The following are global scalars */

#include "../include/macros.h"
#include "../include/su3.h"

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

EXTERN	Real beta;
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  int volume;		/* volume of lattice = nx*ny*nz*nt */
/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */
#define N_POINTERS 8	/* Number of generic pointers */
EXTERN char ** gen_pt[N_POINTERS];

struct site {
	short x,y,z,t;
	char parity;
	int index;
	su3_matrix link[4];
};
typedef struct site site;

EXTERN site *lattice;

#endif /* _LATTICE_H */
