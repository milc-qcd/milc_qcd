#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD QCD program, version 2
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/generic_ks.h" /* For ferm_links_t and ks_action_paths */
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/random.h"    /* For double_prn */
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
	su3_matrix link_save[4];

	/* The Kogut-Susskind phases, which have been absorbed into
	   the matrices.  Also the antiperiodic boundary conditions.  */
	Real phase[4];

	/* 3 element complex vectors */
	su3_vector phi;		/* Gaussian random source vector */
	su3_vector resid;	/* conjugate gradient residual vector */
	su3_vector cg_p;	/* conjugate gradient change vector */
	su3_vector xxx;		/* solution vector = Kinverse * phi */
	su3_vector ttt;		/* temporary vector, for K*ppp */

	su3_vector g_rand[MAX_SRC];	/* Gaussian random vector */
	su3_vector qprop[MAX_SRC];	/* light quark propagators */

	/* temporary vectors and matrices */
	su3_vector tempvec[4];	/* One for each direction */ 
	su3_matrix s_link,s_link_f,t_link_f,diag,staple;
	su3_vector_src resid_src;
	dble_su3_vec_src dtmpvecs[2];


#ifdef RAN_GAUGE
	su3_matrix rgt;
#endif

	/* Additional members of site structure as needed
	   by projects */
#ifdef ADDSITE
#include ADDSITE
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
EXTERN  int volume;			/* volume of lattice = nx*ny*nz*nt */
EXTERN	int iseed;		/* random number seed */
#define MAX_LEVEL 5	/* maximal number of smearing levels */
EXTERN	int no_smear_level, smear_num[MAX_LEVEL], off_axis_flag;
EXTERN  int tot_smear;  /* running total of smearing steps for lattice */
EXTERN	int niter, nrestart, num_src, r0;
EXTERN	Real smear_fac;
EXTERN	Real mass, rsqmin;
EXTERN	char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN	int saveflag;	/* do with lattice: 1=save; */
EXTERN	int total_iters;
EXTERN	int phases_in;	/* 1 if KS phases absorbed into matrices */

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
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8
EXTERN char ** gen_pt[N_POINTERS];

/* Storage for definition of the quark action */
EXTERN ferm_links_t        fn_links;
EXTERN ks_action_paths ks_act_paths;

#endif /* _LATTICE_H */
