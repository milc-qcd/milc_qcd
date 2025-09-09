#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice.

   Directory for dynamical improved KS action.  Allow:
	arbitrary paths in quark action 
	arbitrary paths in gauge action (eg Symanzik imp.)

   If "FN" is defined,
     Includes storage for Naik improvement (longlink[4], templongvec[4],
     gen_pt[16], etc.
     Includes storage for "fat links"  (fatlink[4])
*/

#include "defines.h"
#include "params.h"
#include "../include/generic_quark_types.h"
#include "../include/generic_ks.h" /* For ferm_links_t and ks_action_paths */
#include "../include/random.h"    /* For double_prn */
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */
#include "../include/fermion_links.h"

/* Begin definition of site structure */

#include "../include/su3.h"
#include "../include/random.h"   /* For double_prn */

/* The lattice is an array of sites.  */
#define MOM_SITE   /* If there is a mom member of the site struct */
typedef struct {
    /* The first part is standard to all programs */
	/* coordinates of this site */
	short x,y,z,t;
	/* is it even or odd? */
	char parity;
	/* my index in the array */
	uint32_t index;
#ifdef SITERAND
	/* The state information for a random number generator */
	double_prn site_prn;
	/* align to double word boundary (kludge for Intel compiler) */
	int space1;
#endif

/* ------------------------------------------------------------ */
/*   Now come the physical fields, program dependent            */
/* ------------------------------------------------------------ */
	/* gauge field */
	su3_matrix link[4] ALIGNMENT;	/* the fundamental field */

	/* antihermitian momentum matrices in each direction */
 	anti_hermitmat mom[4] ALIGNMENT;
#ifdef FERMION_FORCE
        su3_matrix ansmom[4];
#endif

	/* The Kogut-Susskind phases, which have been absorbed into 
		the matrices.  Also the antiperiodic boundary conditions.  */
 	Real phase[4];
 	su3_vector g_rand;	/* Gaussian random vector*/
	/* Use trick of combining xxx=D^adj D)^(-1) on even sites with
	   Dslash times this on odd sites when computing fermion force */
} site;

/* End definition of site structure */

/* Definition of globals */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

/* The following are global scalars 
   beta is overall gauge coupling factor
   mass1 and mass2 are quark masses
   u0 is tadpole improvement factor, perhaps (plaq/3)^(1/4)
*/
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  size_t volume;		/* volume of lattice = nx*ny*nz*nt */
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
EXTERN	uint32_t iseed;		/* random number seed */
EXTERN  int nmass;
EXTERN  Real beta;
EXTERN  Real mass,u0;
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, RELOAD_BINARY,
			   RELOAD_CHECKPOINT, FRESH */
EXTERN	int saveflag;	/* do with lattice: FORGET, SAVE, SAVE_BINARY,
			   SAVE_CHECKPOINT */
EXTERN  int savelongflag, savefatflag;
EXTERN  int srcflag, ansflag;
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN  char savelongfile[MAXFILENAME],savefatfile[MAXFILENAME];
EXTERN  char stringLFNlong[MAXFILENAME], stringLFNfat[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN  char srcfile[MAXFILENAME],ansfile[MAXFILENAME];
EXTERN  int inverttype;
EXTERN  int niter;
EXTERN  int nrestart;
EXTERN  Real rsqprop;
EXTERN	int total_iters;
EXTERN	int hisq_svd_counter;
EXTERN	int hisq_force_filter_counter;
EXTERN  int phases_in; /* 1 if KS and BC phases absorbed into matrices */

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	size_t sites_on_node;		/* number of sites on this node */
EXTERN	size_t even_sites_on_node;	/* number of even sites on this node */
EXTERN	size_t odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

EXTERN gauge_file *startlat_p;

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;


/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

EXTERN su3_matrix *ape_links;
EXTERN int refresh_ape_links;
EXTERN int ape_links_ks_phases;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 16
EXTERN char ** gen_pt[N_POINTERS];

/* Storage for definition of the quark action */
EXTERN fermion_links_t        *fn_links;

/* Naik terms */
EXTERN int n_order_naik_total;

/* For eigenpair calculation */
EXTERN int Nvecs_tot;
EXTERN Real *eigVal; /* eigenvalues of D^dag D */
EXTERN su3_vector **eigVec; /* eigenvectors */

EXTERN int n_pseudo;
EXTERN int max_rat_order;

EXTERN int n_order_naik_total;
EXTERN int n_pseudo_naik[MAX_N_PSEUDO];
EXTERN int n_orders_naik[MAX_N_PSEUDO];
#endif /* _LATTICE_H */
