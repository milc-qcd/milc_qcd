#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice_cl.h ********************************/

/* include file for MIMD version 6
   This file defines global scalars and the fields in the lattice. */

/* Modifications:
   8/13/97 wilson_propagator and spin_wilson_vector to su3.h C.D.
   8/10/96  Changed handling of propagator file names; changed same for param C.D.
*/

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

	/* wilson complex vectors */
 	wilson_vector psi;	/* solution vector */
 	wilson_vector chi;	/* source vector */
 /* 	wilson_vector p; conjugate gradient change vector:
				overwrites half the source */
 	wilson_vector mp;	/* another CG vector */
	wilson_vector tmp;	/* another temporary CG vector */
 /*	wilson_vector r; residue: overwrites half the source */
	/* wilson half vector (temporary used in dslash_w_site) */
	half_wilson_vector htmp[MAXHTMP];

	/* storage for one quark_propagator, for four source spins, three source colors */
	wilson_propagator quark_propagator;

	/* Storage for 1/3 quark_propagator, for "rotation" */
	spin_wilson_vector rot_propagator;

  /* extra site members for bicongrad and nondegenerate (HL) cases */
#ifdef BI

#ifdef HL
#include "addsite_clhl_bi.h"
#else
#include "addsite_clov_bi.h"
#endif

#else

#ifdef HL
#include "addsite_clhl_cg.h"
#endif

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
EXTERN	int niter,nrestart,wallflag;
#define MAX_KAP 6
EXTERN	Real kappa,source_r0,kap[MAX_KAP],resid[MAX_KAP];
EXTERN	Real clov_c,u0;
EXTERN	int num_kap;		/* max number of kappa's <= MAX_KAP */
EXTERN	char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN  char savefile_w[MAX_KAP][MAXFILENAME], 
  startfile_w[MAX_KAP][MAXFILENAME];
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

#endif /* _LATTICE_H */
