/****************************** lattice_k.h ********************************/

/* include file for the MIMD heavy-light, code to sum up hopping
	 expansion performed elsewhere  (version 4)
   This file defines global scalars  */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif
#define PI 3.14159265358979323846
#define VERSION_NUMBER 59354  /* For our own binary meson propagator file */
#define WRITELAST 0
#define WRITEALL 1
#define NCHANNELS 4  /* number of meson channels */
#define EPS 1.0e-6


/* The following are global scalars */
EXTERN  int nx,ny,nz,nt;  /* lattice dimensions */
EXTERN  int source_t;
EXTERN	int nhop,writeflag;
EXTERN	Real beta,kappa,kappa_c,kappa_h;
EXTERN	char savefile[73],savefile_k[MAXFILENAME],startfile_m[MAXFILENAME];
EXTERN  int nchannels;
EXTERN	int startflag_m; /* beginning meson hop results: START_ASCII, START_BINARY */
EXTERN	int saveflag_k;	/* save summed meson props: SAVE_ASCII, SAVE_BINARY */

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN  int sites_on_node;    /* number of sites on this node */
EXTERN  int  even_sites_on_node; /* number of even sites on this node */
EXTERN  int  odd_sites_on_node;  /* number of odd sites on this node */
EXTERN  int  number_of_nodes;  /* number of nodes in use */
EXTERN  int  this_node;		/* node number of this node */

/* there is no structure for passing simulation parameters to each node */
/* since only one node does anything */

