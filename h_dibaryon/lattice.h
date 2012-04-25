#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/generic_wilson.h" 
#include "../include/generic_quark_types.h"
#include "../include/random.h"
#include "../include/macros.h"
#include "../include/io_lat.h"    /* For gauge_file */

/* Begin definition of site structure */

#include "../include/su3.h"
#include "../include/complex.h"

typedef struct { su3_vector p[2]; } pauli_vector;
typedef struct { pauli_vector p[2]; } spin_pauli_vector;
typedef struct { spin_pauli_vector c[3]; } pauli_propagator;

typedef struct { complex q[15]; } dipauli_vector;
typedef struct { dipauli_vector q[15]; } dipauli_propagator;


/* The lattice is an array of sites.  Within each node the even sites will
   be stored first, then the odd sites.
*/


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

	/* wilson complex vectors */
 	wilson_vector psi;	/* solution vector */
 	wilson_vector chi;	/* source vector */
 /* 	wilson_vector p; conjugate gradient change vector:
				overwrites half the source */
 	wilson_vector mp;	/* another CG vector */
	wilson_vector tmp;	/* another temporary CG vector */
 /*	wilson_vector r; residue: overwrites half the source */
	/* wilson half vector (temporary used in dslash_w_site) */
	half_wilson_vector htmp[2];

	/* storage for one quark_propagator, for four source spins, three source colors */
	wilson_propagator quark_propagator;
	wilson_propagator quark_prop2;
	/* Two nonrelativistic propagators */
	pauli_propagator nr_prop1;
	pauli_propagator nr_prop2;

	/* Two diquark propagators */
	dipauli_propagator diquark_prop1;
	dipauli_propagator diquark_prop2;

#ifdef BI
#include "addsite_clhl_bi.h"
#endif

} site;

/* End definition of site structure */

/* Definition of globals */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

#define GAUGE_FIX_TOL 1.0e-7 /* For gauge fixing */

/* The following are global scalars */
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  int volume;		/* volume of lattice = nx*ny*nz*nt */
EXTERN	int niter,nrestart,wallflag;
EXTERN	int source_time;       /* timeslice of the initial source */
EXTERN	Real kappa,kappa_heavy,kappa_light;
#define MAX_KAP 6
EXTERN  Real width,source_r0,kap[MAX_KAP],resid[MAX_KAP];
EXTERN	Real clov_c,u0;
EXTERN  int num_kap;            /* total number of kappas <= MAX_KAP */
EXTERN	int num_kap_heavy;	/* number of heavy kappas <= MAX_KAP */
EXTERN	int num_kap_light;	/* number of light kappas <= MAX_KAP */
EXTERN  int nr_forw_back;       /* NR source forw (1) backw (2) both (3) */
EXTERN	char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN  char savefile_w[MAX_KAP][MAXFILENAME], startfile_w[MAX_KAP][MAXFILENAME];
EXTERN  char scratchstem_w[MAXFILENAME];
EXTERN  int scratchflag;        /* scratch file mode: SAVE_SERIAL, SAVE_CHECKPOINT */
EXTERN	int startflag;		/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN  int fixflag;  /* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */
EXTERN	int saveflag;		/* save lattice: SAVE_ASCII, SAVE_BINARY */
EXTERN	int startflag_w[MAX_KAP];  /* beginning wilson: 
			   RELOAD_ASCII, RELOAD_BINARY, RELOAD_PARALLEL */
EXTERN	int saveflag_w[MAX_KAP];   /* save propagator: SAVE_ASCII, SAVE_BINARY
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

EXTERN quark_source wqs[MAX_KAP];
EXTERN quark_source wqstmp;  /* Temporary */

EXTERN quark_invert_control qic;
EXTERN dirac_clover_param dcp;

EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;


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
