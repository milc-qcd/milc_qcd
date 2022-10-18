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
EXTERN  size_t volume;		/* volume of lattice = nx*ny*nz*nt */
/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	size_t sites_on_node;		/* number of sites on this node */
EXTERN	size_t even_sites_on_node;	/* number of even sites on this node */
EXTERN	size_t odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
#define N_POINTERS 8	/* Number of generic pointers */
EXTERN char ** gen_pt[N_POINTERS];
/* Further source description */
#ifdef CONTROL
int n_spins = 4;                /* Number of spins generated */
int source_loc[4] = {0,0,0,0};  /* Source location */
int spins[4] = {0,1,2,3};       /* List of spins generated */
#else
extern int n_spins;
extern int source_loc[4];
extern int spins[4];
#endif

struct site {
	short x,y,z,t;
	char parity;
	uint32_t index;
	double_prn site_prn;
#ifndef NO_GAUGE_FIELD
	su3_matrix link[4] ALIGNMENT;
#endif
};
typedef struct site site;

EXTERN Real boundary_phase[4];
EXTERN site *lattice;

EXTERN su3_matrix *ape_links;
EXTERN int refresh_ape_links;

#endif /* _LATTICE_H */
