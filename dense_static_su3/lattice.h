#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD QCD program, version 4
   This file defines global scalars and the fields in the lattice. */


#include "defines.h"
#include "../include/random.h"
#include "../include/macros.h"
#include "../include/io_lat.h"    /* For gauge_file */

/* Begin definition of site structure */

#include "../include/su3.h"
#include "../include/complex.h"
#include "../include/random.h"   /* For double_prn */

/* =============================================================*/
/* --------------------The lattice -----------------------------*/
/* =============================================================*/
typedef struct {
/* --------------------The first part is standard to all programs */
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

/* ---------------------The physical fields, program dependent */
	/* gauge field */
	su3_matrix link[4];  /* connection */

	/* temporary matrices */
	/* P-loop */
	complex ploop;  /* holds tr P */
        su3_matrix ploop_t;  /* holds PU^*, "P-loop - link" */ 

	/* workspace matrices */
	su3_matrix staple;
	su3_matrix tempmat1;
	su3_matrix tempmat2;

	complex qprob[4];	/* probability of N quarks at this site */

#ifdef LIGHT_PBP	/* for light quark psi-bar-psi measurement */
        /* The Kogut-Susskind phases, which have been absorbed into
                the matrices.  Also the antiperiodic boundary conditions.  */
        Real phase[4];
        /* 3 element complex vectors */
        su3_vector phi;         /* Gaussian random source vector */
        su3_vector resid;       /* conjugate gradient residual vector */
        su3_vector cg_p;        /* conjugate gradient change vector */
        su3_vector xxx;         /* solution vector = Kinverse * phi */
        su3_vector ttt;         /* temporary vector, for K*ppp */
        su3_vector g_rand;      /* Gaussian random vector*/
        su3_vector tempvec[4];  /* One for each direction */
#endif /* LIGHT_PBP */
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
EXTERN	int warms,trajecs,steps_over,steps_update,propinterval;
EXTERN	Real beta;
EXTERN	Real C;
EXTERN	Real epsilon;
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN  int fixflag;  /* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */
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

#ifdef LIGHT_PBP
/* a few extra variables for light quark psi-bar-psi measurements */
EXTERN  Real rsqmin,rsqprop,mass;
EXTERN  int nflavors,total_iters;
#endif /* LIGHT_PBP */

EXTERN gauge_file *startlat_p;

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;






/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

/*===============================================================*/
/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8
EXTERN char ** gen_pt[N_POINTERS];

/* Some general purpose macros */
/* under construction
#define DEBUG_ON #ifdef DEBUG
#define DEBUG_OFF #endif
*/

/* Storage for definition of the quark action */
EXTERN ferm_links_t        fn_links;
EXTERN ks_action_paths ks_act_paths;

#endif /* _LATTICE_H */

