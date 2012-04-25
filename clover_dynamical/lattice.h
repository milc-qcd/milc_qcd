#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD QCD program, version 2
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */
#include "../include/random.h"    /* For double_prn */
#include "../include/generic_clover.h" /* For clover */

/* Begin definition of site structure */

#include "../include/su3.h"
#include "../include/random.h"   /* For double_prn */

/* The lattice is an array of sites.  */
#define MOM_SITE   /* If there is a mom member of the site struct */
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
#ifdef HMC_ALGORITHM
 	su3_matrix old_link[4];
	/* For accept/reject */
#endif

	/* antihermitian momentum matrices in each direction */
 	anti_hermitmat mom[4];

	/* Wilson complex vectors */
	wilson_vector g_rand;	/* gaussian random vector */
	wilson_vector psi;	/* solution vector */
        wilson_vector chi;	/* source vector */
#ifdef BI
        wilson_vector sss;      /* internal biconjugate gradient vector */
        wilson_vector ttt;      /* internal biconjugate gradient vector */
        wilson_vector rv;       /* internal biconjugate gradient vector */
        wilson_vector v;        /* internal biconjugate gradient vector */
#endif
        wilson_vector p;        /* conjugate gradient change vector */
        wilson_vector mp;       /* another CG vector */
        wilson_vector tmp;      /* another temporary CG vector */
	wilson_vector r; 	/* residue */
        /* wilson half vector (temporary used in dslash_w_site) */
        half_wilson_vector htmp[8];
/**half_wilson_vector htmp2[8];**/ /* TEMP FOR TESTING */
#ifdef PHI_ALGORITHM
 	wilson_vector old_psi;	/* For predicting next psi */
#endif
#ifdef SPECTRUM
	wilson_matrix quark_propagator;
	    /* For four source spins, three source colors */
	wilson_matrix rotated_propagator;
	    /* For clover-rotated operators */
#endif
	/* temporary vectors and matrices */
	su3_matrix tempmat1,tempmat2,staple;

    /* align to double word boundary */
	/**double space2;**/
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
EXTERN	int warms,trajecs,steps,niter,nrestart,propinterval,nflavors;
EXTERN	Real rsqmin,rsqprop,beta,kappa,clov_c,u0;
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

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */


/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;

#define nloop 3
#define nreps 1
#ifdef GB
#define nloop_m 4
#define nreps_m 3
#endif
#define max_num 400

/* global defns for general action (check this works...) */
EXTERN int loop_ind[nloop][10], loop_length[nloop];
EXTERN int loop_table[nloop][max_num][10],loop_num[nloop],loop_char[max_num];
EXTERN Real loop_coeff[nloop][nreps];
EXTERN int loop_ch[nloop][max_num],ch;

EXTERN Real loop_term[48][nreps];

#ifdef GB
EXTERN int loop_ind_m[nloop_m][10], loop_length_m[nloop_m];
EXTERN int loop_table_m[nloop_m][max_num][10],loop_num_m[nloop_m],
loop_char_m[max_num];
EXTERN int loop_ch_m[nloop_m][max_num];
#endif

/* definitions for arccos */
#define NMAX 401
 
EXTERN Real acos_table[NMAX];
EXTERN Real acos_deriv[NMAX];
EXTERN Real delta_acos;

EXTERN quark_invert_control qic;
EXTERN dirac_clover_param dcp;

EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8	/* Number of generic pointers */
/* NEED 8 WHEN GAUGEFIXING */
EXTERN char ** gen_pt[N_POINTERS];

/* Storage for the clover term */
EXTERN clover *gen_clov;

#endif /* _LATTICE_H */
