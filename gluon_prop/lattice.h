#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/generic_ks.h" /* For ferm_links_t and ks_action_paths */
#include "../include/macros.h"  /* For MAXFILENAME */
#include "../include/io_lat.h"	/* For gauge_file */
#ifdef QUARK_PROP
#include "../include/fermion_links.h" /* For fermion_links_t */
#endif

/* Begin definition of site structure */

#include "../include/su3.h"
#include "../include/random.h"	/* For double_prn */

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

#ifdef GLUON_PROP
	/* vector potential */
	su3_matrix a_mu[4];
	/* propagator */
	Real scalar_prop, long_prop;
#endif

#ifdef QUARK_PROP
	/* propagator */
	complex trace_prop[MAX_NUM_MASS];
	Real qprop[2];
	Real sum_qprop[2];

	/* The Kogut-Susskind phases, which have been absorbed into 
		the matrices.  Also the antiperiodic boundary conditions.  */
	Real phase[4];

	/* 3 element complex vectors */
	su3_vector phi;		/* Gaussian random source vector */
	su3_vector resid;	/* conjugate gradient residual vector */
	su3_vector cg_p;	/* conjugate gradient change vector */
	su3_vector xxx1;	/* solution vector = Kinverse * phi */
	su3_vector ttt;		/* temporary vector, for K*ppp */

	/* temporary vectors and matrices */
	su3_vector tempvec[4];	/* One for each direction */
#ifdef FN
	su3_vector templongvec[4];	/* One for each direction */
	su3_vector templongv1;
#endif
	/* temporary matrix (Note: staple is used in quark_stuff.c) */
	su3_matrix staple;
#endif

        /* temporary matrix   (Note: mom is used in map_milc_to_qop.c) */
        /* We really need to get rid of these */
        su3_matrix mom[4];
	/* temporary matrices */
	su3_matrix tempmat1, tempmat2;
#ifdef GLUON_PROP
#ifdef GFIX
	su3_vector tempvec[1];
	su3_matrix staple;
#endif
#endif

	/* temporary Real */
	Real tempfloat;

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
EXTERN	int iseed;		/* random number seed */
EXTERN	int niter;		/* max number of CG iterations */
EXTERN	int nrestart;		/* max number of CG restarts */
EXTERN	int num_mass;		/* number of masses */
EXTERN	Real beta,quarkmass,u0;
EXTERN	Real mass[MAX_NUM_MASS];
EXTERN	Real rsqprop;
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN	char stringLFN[MAXFILENAME];	/** ILDG LFN if applicable **/
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, RELOAD_BINARY,
			   RELOAD_CHECKPOINT, FRESH */
EXTERN  int fixflag;	/* gauge fix: COULOMB_GAUGE_FIX, LANDAU_GAUGE_FIX,
			   NO_GAUGE_FIX */
EXTERN  int fixflag_ft;	/* FFT gauge fix: same as above */
EXTERN	int saveflag;	/* do with lattice: FORGET, SAVE, SAVE_BINARY,
			   SAVE_CHECKPOINT */
EXTERN	char ksstartfile[MAX_NUM_MASS][MAXFILENAME];	/* from where to read KS prop */
EXTERN	char kssavefile[MAX_NUM_MASS][MAXFILENAME];	/* where to store KS prop */
EXTERN	int ksstartflag[MAX_NUM_MASS];	/* read KS propagator or not:
			   CONTINUE, FRESH, RELOAD_ASCII, RELOAD_SERIAL */
EXTERN	int kssaveflag[MAX_NUM_MASS];	/* save KS propagator or not:
			   FORGET, SAVE_ASCII, SAVE_SERIAL, SAVE_SERIAL_FM,
			   SAVE_SERIAL_TSLICE */
EXTERN	int run_CG_flag[MAX_NUM_MASS];	/* Do the inversion, or not? */
EXTERN	int total_iters;
EXTERN	int phases_in;	/* 1 if KS and BC phases absorbed into matrices */

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN	int this_node;		/* node number of this node */

EXTERN	gauge_file *startlat_p;


/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN  double_prn node_prn ;

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

EXTERN su3_matrix *ape_links;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 16
EXTERN char ** gen_pt[N_POINTERS];

#ifdef QUARK_PROP
/* Storage for definition of the quark action */
EXTERN fermion_links_t *fn_links;
#endif


#endif /* _LATTICE_H */
