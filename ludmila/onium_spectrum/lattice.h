#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice_cl.h ********************************/

/* include file for MIMD version 6
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/random.h"   /* For double_prn */
#include "../include/macros.h"   /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */

/* Begin definition of site structure */

#include "../include/su3.h"

typedef struct {
  /* The first part is standard to all programs */
	/* coordinates of this site */
  short x,y,z,t;
  /* is it even or odd? */
  char parity;
  /* my index in the array */
  int index;
  
  /* Now come the physical fields, program dependent */
  /* gauge field */
  su3_matrix link[4];
  
  /* wilson complex vectors - removed*/
	wilson_vector psi;	/* solution vector */
 	wilson_vector chi;	/* source vector */

  /* 	wilson_vector p; conjugate gradient change vector:
	overwrites half the source */

 	wilson_vector mp;	/* another CG vector */
	wilson_vector tmp;	/* another temporary CG vector */
 /*	wilson_vector r; residue: overwrites half the source */ 
  /* wilson half vector (temporary used in dslash_ks???) */
  /*half_wilson_vector htmp[MAXHTMP];*/
  
  /* storage for one quark_propagator, for four source spins, three source colors */
  wilson_propagator quark_propagator;
  wilson_propagator quark_propagator_copy;
  wilson_propagator antiquark_propagator;
  wilson_propagator antiquark_propagator_copy;  
  /* Storage for 1/3 quark_propagator, for "rotation" */
  spin_wilson_vector rot_propagator;
  
  half_wilson_vector htmp[MAXHTMP];	
  
  /* storage for staggered propagator*/
  
  su3_matrix stag_propagator;
  su3_matrix stag_propagator_copy;  
  /* smearing function */
  
  complex w,w1,w2;

  /* The Kogut-Susskind phases, which have been absorbed into 
     the matrices.  Also the antiperiodic boundary conditions.  */
  float phase[4];
  
  /* 3 element complex vectors */
  su3_vector phi;		/* Gaussian random source vector */
  su3_vector resid;	/* conjugate gradient residual vector */
  su3_vector cg_p;	/* conjugate gradient change vector */
  su3_vector xxx;		/* solution vector = Kinverse * phi */
  su3_vector ttt;		/* temporary vector, for K*ppp */
  su3_vector g_rand;	/* Gaussian random vector*/
#ifdef SITERAND
  	/* The state information for a random number generator */
  	double_prn site_prn;
		/* align to double word boundary (kludge for Intel compiler) */
		int space1;
#endif
  /* temporary vectors and matrices */
  su3_vector tempvec[4];	/* One for each direction */
  su3_matrix tempmat1;
  
  
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
EXTERN	int niter,nrestart,wallflag;
#define MAX_KAP 6
EXTERN	float kappa,source_r0,kap[MAX_KAP],resid[MAX_KAP];
EXTERN	float clov_c,u0;
EXTERN	int num_kap;		/* max number of kappa's <= MAX_KAP */
EXTERN	char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN  char savefile_w[MAX_KAP][MAXFILENAME], 
  startfile_w[MAX_KAP][MAXFILENAME], a_startfile_w[MAXFILENAME];
EXTERN  char start_ks_prop_file[MAXFILENAME];
EXTERN  char scratchstem_w[MAXFILENAME];
EXTERN  int scratchflag;        /* scratch file mode: SAVE_SERIAL, SAVE_CHECKPOINT */
EXTERN	int startflag;		/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN  int fixflag;  /* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */
EXTERN	int saveflag;		/* save lattice: SAVE_ASCII, SAVE_BINARY */
EXTERN	int startflag_w[MAX_KAP];  /* beginning wilson: 
			   RELOAD_ASCII, RELOAD_BINARY, RELOAD_PARALLEL */
EXTERN	int saveflag_w[MAX_KAP];   /* save propagator: SAVE_ASCII, SAVE_BINARY*/

EXTERN  char smearfile[MAX_KAP][MAXFILENAME];
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

EXTERN wilson_quark_source wqs[MAX_KAP];
EXTERN wilson_quark_source wqstmp;  /* Temporary */

EXTERN quark_invert_control qic;
EXTERN dirac_clover_param dcp;

EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN  double_prn node_prn ;

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN site *lattice;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8	/* Number of generic pointers */
/* NEED 8 WHEN GAUGEFIXING */
EXTERN char ** gen_pt[N_POINTERS];
EXTERN	int iseed;

EXTERN  float d1[MAX_KAP];

EXTERN  int num_smear;
EXTERN  int format[MAX_KAP];
EXTERN  int a_format;
EXTERN	float rsqmin,rsqprop,beta,mass;
EXTERN	int warms,trajecs,steps,niter,propinterval,nflavors;
EXTERN	float epsilon;
EXTERN  int phases_in; /* 1 if KS and BC phases absorbed into matrices */
EXTERN  float a_d1;
#endif /* _LATTICE_H */
