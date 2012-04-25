#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */
#include "../include/generic_wilson.h" 

/* does NOT use generic_clover routines! */
typedef struct { complex tr[2][15]; } triangular;
typedef struct { Real di[2][6]; } diagonal;

/* Begin definition of site structure */

#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/random.h"   /* For double_prn */
#include "../include/generic_quark_types.h"

typedef struct {
    /* The first part is standard to all programs */
	/* coordinates of this site */
	short x,y,z,t;
	/* is it even or odd? */
	char parity;
	/* my index in the array */
	int index;
#ifdef SITERAND
        /* The state information for a random number generator */
        double_prn site_prn;
        /* align to double word boundary (kludge for Intel compiler) */
        int space1;
#endif

    /* Now come the physical fields, program dependent */
	/* gauge field */
	su3_matrix link[4];

        su3_matrix blocked_link[NLINK];

        su3_matrix tempmat1,tempmat2;

	/* wilson complex vectors */
 	wilson_vector psi;	/* solution vector */
 	wilson_vector chi;	/* source vector */
  	wilson_vector p; /* conjugate gradient change vector */
 	wilson_vector mp;	/* another CG vector */
 	wilson_vector r; /* residue:  */
	/* wilson vector (temporary used in delta0) */
	wilson_vector htmp[2];
/* Bi-CG code */
        wilson_vector sss;      /* internal biCG vector */
        wilson_vector ttt;      /* internal biCG vector */
        wilson_vector rv;       /* internal biCG vector */
        wilson_vector v;        /* internal biCG vector */
/* clover stuff */
        wilson_vector tmp;  /* another temporary CG vector used in clover */
        su3_matrix staple;
        triangular clov;
        diagonal clov_diag;
#ifdef RANDOM
        wilson_vector source;      /* to save the random source */
#endif
#ifdef CVC
complex cvc[4],cvct[4];
complex hvc[4],hvct[4];
#endif
#ifdef SPECTRUM
	/* storage for one quark_propagator, for four source spins, three source colors */
	wilson_propagator quark_propagator;
        /* Storage for 1/3 quark_propagator, for "rotation" etc*/
	spin_wilson_vector extra_propagator;
#endif
} site;

/* End definition of site structure */

/* Definition of globals */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

/* The following are global scalars */
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  int volume;		/* volume of lattice = nx*ny*nz*nt */
#ifdef RANDOM
EXTERN  int iseed;              /* random number seed */
#endif
EXTERN	int niter,nrestart,wallflag;
#define MAX_MASSES 6
EXTERN	Real source_r0,mass[MAX_MASSES],resid[MAX_MASSES];
EXTERN	int num_masses;		/* max number of masses <= MAX_MASSES */
EXTERN  Real m0,lambda[5],rho[5];
EXTERN  char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN  char savefile_w[MAX_MASSES][MAXFILENAME],
             startfile_w[MAX_MASSES][MAXFILENAME];
EXTERN  char scratchstem_w[MAXFILENAME];
EXTERN  int scratchflag;/* scratch file mode: SAVE_SERIAL, SAVE_CHECKPOINT */


EXTERN	int startflag; 	/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN  int fixflag;  /* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */
EXTERN	int saveflag;		/* save lattice: SAVE_ASCII, SAVE_BINARY */
EXTERN	int startflag_w[MAX_MASSES];  /* beginning wilson: 
			   RELOAD_ASCII, RELOAD_BINARY, RELOAD_PARALLEL */
EXTERN	int saveflag_w[MAX_MASSES]; /* save propagator: SAVE_ASCII, SAVE_BINARY
				   SAVE_PARALLEL */
EXTERN	int total_iters;

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

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */


EXTERN quark_source wqs[MAX_MASSES];
EXTERN quark_source wqstmp;  /* Temporary */




EXTERN Real clover_term;

/* information for derivative and link terms */
/*number of paths in the hypercube for 1 and gamma-mu paths */
#define NLINK 40
EXTERN int offset[NLINK][4];
EXTERN int label[NLINK],off_max;



/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN  double_prn node_prn ;

/* The lattice is an array of sites.  Within each node the even sites will
   be stored first, then the odd sites.
*/


/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8	/* Number of generic pointers */
/* NEED 8 WHEN GAUGEFIXING */
EXTERN char ** gen_pt[N_POINTERS];
#endif /* _LATTICE_H */


