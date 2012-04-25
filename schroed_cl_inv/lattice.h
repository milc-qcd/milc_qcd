#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/generic_wilson.h" 
#include "../include/generic_quark_types.h"
#include "../include/macros.h"	/* For MAXFILENAME, EVENFIRST */
#include "../include/random.h"	/* For double_prn */
#include "../include/io_lat.h"	/* For gauge_file */
#include "../include/generic_clover.h" /* For clover */

/* Begin definition of site structure */

#include "../include/su3.h"
#include "../include/random.h"	/* For double_prn */

typedef struct {
  int color;
  int spin;
  int type;
  double kappa;
  char descrp[128];
} my_quark_source;

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
	su3_matrix boundary[3];	/* spatial boundary links (for t=nt) */
	su3_matrix link_tmp[4];

	/* wilson complex vectors */
 	wilson_vector psi;	/* solution vector */
 	wilson_vector chi;	/* source vector */
#ifdef BI
	wilson_vector sss;	/* internal biconjugate gradient vector */
#endif
 	wilson_vector mp;	/* another CG vector */
	wilson_vector tmp;	/* another temporary CG vector */
	wilson_vector tmpb;	/* auxiliary */
	/* wilson half vector (temporary used in dslash_w_site) */
	half_wilson_vector htmp[2];

	/* storage for one full quark_propagator, for forward propagators */
	wilson_propagator forw_quark_prop;

	/* storage for one full quark_propagator, for backward propagators */
	wilson_propagator backw_quark_prop;

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
EXTERN	int niter,nrestart;
#define MAX_KAP 6
EXTERN	Real kappa,kap[MAX_KAP],resid[MAX_KAP];
EXTERN	Real clov_c;
EXTERN	int num_kap;		/* max number of kappa's <= MAX_KAP */
EXTERN	int bc_flag;		/* flag for gauge field bc */
EXTERN	Real ferm_phases[3];	/* fermion phase factors */
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN	int startflag;		/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN	int total_iters;
EXTERN	int num_smear;
EXTERN	Real alpha;		/* APE smearing parameter (Boulder convention) */
EXTERN	Real beta;
EXTERN	Real c_t11;

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */
EXTERN  my_quark_source wqs[MAX_KAP];
EXTERN	quark_invert_control qic;
EXTERN	dirac_clover_param dcp;

EXTERN	gauge_file *startlat_p;

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN	double_prn node_prn ;

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN	site *lattice;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8	/* Number of generic pointers */
/* NEED 8 WHEN GAUGEFIXING */
EXTERN	char ** gen_pt[N_POINTERS];

/* Storage for the clover term */
EXTERN clover *gen_clov;

#endif	/* _LATTICE_H */
