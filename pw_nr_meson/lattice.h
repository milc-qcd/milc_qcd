#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice_cl.h ********************************/

/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/generic_wilson.h"
#include "../include/random.h"   /* For double_prn */
#include "../include/macros.h"   /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */
#include "pauli_prop.h"
#include "../include/generic_clover.h" /* For clover */

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

#ifdef SITERAND
  /* The state information for a random number generator */
  double_prn site_prn;
  /* align to double word boundary (kludge for Intel compiler) */
  int space1;
#endif

  half_wilson_vector htmp[MAXHTMP];
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
EXTERN	Real kappa;
EXTERN	Real clov_c,u0;
EXTERN	Real resid,relresid;
EXTERN	int niter, nrestart;
EXTERN  int fixflag;
EXTERN	char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char startfile_w[NSM][MAXDIR][MAXFILENAME], a_startfile_w[MAXFILENAME];
EXTERN  char savefile_w[NSM][MAXDIR][MAXFILENAME], a_savefile_w[MAXFILENAME];
EXTERN  wilson_quark_source source_wqs[NSM];
EXTERN  wilson_quark_source sink_wqs[NSM];
EXTERN  wilson_quark_source a_wqs;
EXTERN  wilson_quark_source wqstmp;  /* Temporary */
EXTERN  int a_startflag_w;
EXTERN  int a_saveflag_w;
EXTERN	int startflag_w[NSM][MAXDIR];  
EXTERN	int saveflag_w[NSM][MAXDIR];   
EXTERN  char scratchstem_w[MAXFILENAME];
EXTERN  int scratchflag;        /* scratch file mode: SAVE_SERIAL, SAVE_CHECKPOINT */
EXTERN	int startflag;		/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN  int fixflag;  /* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */
EXTERN	int saveflag;		/* save lattice: SAVE_ASCII, SAVE_BINARY */

EXTERN  char smearfile[NSM][MAXFILENAME];
EXTERN  char source_wf_label[NSM][MAXFILENAME];
EXTERN  char sink_wf_label[NSM][MAXFILENAME];

EXTERN  char a0_file[MAXFILENAME];
EXTERN  char b1_file[MAXFILENAME];
EXTERN  char a1_file[MAXFILENAME];
EXTERN  char a2_t2_file[MAXFILENAME];
EXTERN  char a2_e_file[MAXFILENAME];
EXTERN  char a2_file[MAXFILENAME];

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
EXTERN  int sequence;

EXTERN quark_invert_control qic;
EXTERN dirac_clover_param dcp;

EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

EXTERN quark_invert_control qic;
EXTERN dirac_clover_param dcp;

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

//EXTERN  Real d1[MAX_KAP];

EXTERN  int num_smear;

EXTERN int flag;

EXTERN  block_pauli_propagator *antiquark_prop;
EXTERN  block_pauli_propagator *quark_prop;
EXTERN  block_pauli_propagator *quark_prop_smear;

/* Storage for the clover term */
EXTERN clover *gen_clov;

#endif /* _LATTICE_H */
