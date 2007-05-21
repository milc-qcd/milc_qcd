/****************************** lattice_k.h ********************************/

/* include file for the MIMD heavy-light, code to sum up hopping
	 expansion performed elsewhere  (version 5)
   This file defines global scalars  */

#include "defines.h"
#include "../../include/macros.h"   /* For MAXFILENAME */

/* Begin definition of site structure */


/*** A set of dummy definition required to compile comdefs.h and some IO routines***/

typedef struct
{
  short x,y,z,t;
  /* is it even or odd? */
  char parity;
  int index ;
} site ;



/* End definition of site structure */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

#define WRITELAST 0
#define WRITEALL 1
#define NCHANNELS 4  /* number of meson channels */
#define EPS 1.0e-6

/* Definition of globals */

/*enum io_options { RELOAD_ASCII = 90 , RELOAD_BINARY  } ; */
#define RELOAD_BINARY 90

/* The following are global scalars */
EXTERN  int nx,ny,nz,nt;  /* lattice dimensions */
EXTERN  int source_t;
EXTERN	int nhop,writeflag;
EXTERN	Real beta,kappa,kappa_c,kappa_h;
EXTERN	char savefile[73],savefile_k[MAXFILENAME],startfile_m[MAXFILENAME];
EXTERN  int nchannels;  /* = NCHANNELS +2 if extra_sink chosen */
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

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */

/* EXTERN site *lattice; */

