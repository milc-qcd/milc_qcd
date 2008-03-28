#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice_cl.h ********************************/

/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/generic_ks.h" 
#include "../include/generic_wilson.h" 
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
  
  /* Storage for 1/3 quark_propagator, for "rotation" */
  spin_wilson_vector rot_propagator;
  
  half_wilson_vector htmp[MAXHTMP];	
  
  /* storage for staggered propagator*/
  
  su3_matrix stag_propagator;
  su3_vector prop[3];
  su3_matrix stag_propagator_copy;  

  /* smearing function */
  
  complex w,w1,w2;

  /* The Kogut-Susskind phases, which have been absorbed into 
     the matrices.  Also the antiperiodic boundary conditions.  */
  double phase[4];
  
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
  
  /* for the baryon code */
  su3_matrix stag_strange_propagator;
  su3_matrix stag_light_propagator;
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
EXTERN  char job_id[MAXFILENAME];
EXTERN  int volume;		/* volume of lattice = nx*ny*nz*nt */
EXTERN	int nrestart,wallflag;
#define MAX_KAP 6
EXTERN	double kappa,source_r0,kap[MAX_KAP],resid[MAX_KAP];
EXTERN  char kap_label[MAX_KAP][32];
EXTERN	double clov_c,u0;
EXTERN	int num_kap;		/* max number of kappa's <= MAX_KAP */
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, RELOAD_BINARY,
			   RELOAD_CHECKPOINT, FRESH */
EXTERN	int saveflag;	/* do with lattice: FORGET, SAVE, SAVE_BINARY,
			   SAVE_CHECKPOINT */
EXTERN  int fixflag;  /* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN  char savefile_w[MAX_KAP][MAXFILENAME];
EXTERN  char startfile_w[MAX_KAP][MAXFILENAME];
EXTERN  char src_label_w[MAX_KAP][16];
EXTERN  char savefile_c[MAXFILENAME];
EXTERN  char start_ks_prop_file[MAXFILENAME];
EXTERN  char sink_label[MAX_KAP][16];
EXTERN  char scratchstem_w[MAXFILENAME];
EXTERN  int scratchflag;        /* scratch file mode: SAVE_SERIAL, SAVE_CHECKPOINT */
EXTERN	int startflag_w[MAX_KAP];  /* beginning wilson: 
			   RELOAD_ASCII, RELOAD_BINARY, RELOAD_PARALLEL */
EXTERN	int saveflag_w[MAX_KAP];   /* save propagator: SAVE_ASCII, SAVE_BINARY*/
EXTERN	int saveflag_c;   /* save propagator: SAVE_ASCII, SAVE_BINARY*/
EXTERN  int ks_prop_startflag;

EXTERN  char smearfile[MAX_KAP][MAXFILENAME];
EXTERN	int total_iters;
EXTERN  double d1[MAX_KAP];  /*rotation parameter*/
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

EXTERN ks_quark_source ksqs;
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

EXTERN  int log_correlators;
EXTERN	double rsqmin,rsqprop,beta,mass;
EXTERN  char mass_label[32];
EXTERN	int warms,trajecs,steps,niter,propinterval,nflavors;
EXTERN	double epsilon;
EXTERN  int phases_in; /* 1 if KS and BC phases absorbed into matrices */
EXTERN  int num_smear;

/* Variables for the baryon code */

#define MAX_STRANGE 4
#define MAX_LIGHT 4

EXTERN  int num_strange;
EXTERN  int num_light;
EXTERN  char start_ks_strange_file[MAX_STRANGE][MAXFILENAME];
EXTERN  char start_ks_light_file[MAX_LIGHT][MAXFILENAME];
EXTERN  int start_ks_strange_flag[MAX_STRANGE];
EXTERN  int start_ks_light_flag[MAX_LIGHT];
EXTERN  Real m_light[MAX_LIGHT];
EXTERN  Real m_strange[MAX_LIGHT];


#endif /* _LATTICE_H */
