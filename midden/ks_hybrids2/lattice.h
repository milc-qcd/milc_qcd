#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD version 6
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/random.h"    /* For double_prn */
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */

/* Begin definition of site structure */

#include "../include/su3.h"
#include "../include/random.h"   /* For double_prn */

/* The lattice is an array of sites.  */
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
	su3_matrix smearlink[4];
  /* For parallel I/O routines until we change them C.D. */

	/* The Kogut-Susskind phases, which have been absorbed into 
		the matrices.  Also the antiperiodic boundary conditions.  */
 	Real phase[4];

	/* 3 element complex vectors */
 	su3_vector phi;		/* Gaussian random source vector */
 	su3_vector resid;	/* conjugate gradient residual vector */
 	su3_vector cg_p;	/* conjugate gradient change vector */
 	su3_vector xxx;		/* solution vector = Kinverse * phi */
 	su3_vector ttt;		/* temporary vector, for K*ppp */
 	su3_vector g_rand;	/* Gaussian random vector*/
	/* temporary vectors and matrices */
	su3_vector tempvec[4];	/* One for each direction */
	su3_matrix tempmat1;	/* for gaugefix,plaq.[23],ploop[23],barcorr */
	su3_matrix tempmat2;	/* for plaq.[23],ploop[23],barcorr,spect... */

	/* components of F_mu_nu, use defines to give  components names */
	su3_matrix field_strength[6];
	/* quark and antiquark propagators */
	su3_vector quark_source;
	su3_vector quark_prop;
	su3_vector anti_prop;

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
EXTERN  int volume;			/* volume of lattice = nx*ny*nz*nt */
EXTERN	int iseed;		/* random number seed */
EXTERN	int niter;
EXTERN	Real rsqprop,beta,mass;
EXTERN  Real staple_weight;	/* for smearing */
EXTERN  int source_start, source_inc, n_sources;
        /* source time, increment for it, and number of source slices */
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN	int saveflag;	/* do with lattice: 1=save; */
EXTERN	int total_iters;
EXTERN  int phases_in; /* 1 if KS and BC phases absorbed into matrices */

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

EXTERN gauge_file *startlat_p;


/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;


/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN site *lattice;

/* defines for index on field_strength */
#define FS_XY 0
#define FS_XZ 1
#define FS_YZ 2
#define FS_XT 3
#define FS_YT 4
#define FS_ZT 5

#define ON 1	/* for tracking whether phases are in */
#define OFF 0

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8
EXTERN char ** gen_pt[N_POINTERS];
EXTERN char ** gen_pt2[N_POINTERS];


#endif /* _LATTICE_H */


