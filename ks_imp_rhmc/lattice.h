#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/
/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice.

   Directory for dynamical improved KS action.  Allow:
	arbitrary paths in quark action 
	arbitrary paths in gauge action (eg Symanzik imp.)

   If "FN" is defined,
     Includes storage for Naik improvement (longlink[4], templongvec[4],
     gen_pt[16], etc.
*/

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/generic_ks.h" /* For ferm_links_t and ks_action_paths */
#include "../include/random.h"    /* For double_prn */
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */
#include "params.h"
#include "../include/fermion_links.h"  /* For fermion_links_t */
#ifdef HAVE_U1
#include "../include/io_u1lat.h"  /* For gauge_file */
#include "../include/complex.h"
#endif

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
	uint32_t index;
#ifdef SITERAND
	/* The state information for a random number generator */
	double_prn site_prn;
	/* align to double word boundary (kludge for Intel compiler) */
	int space1;
#endif

/* ------------------------------------------------------------ */
/*   Now come the physical fields, program dependent            */
/* ------------------------------------------------------------ */

        /* gauge field */
        su3_matrix link[4] ALIGNMENT;	/* the fundamental field */

	/* antihermitian momentum matrices in each direction */
        anti_hermitmat mom[4] ALIGNMENT;

#ifdef HMC
 	su3_matrix old_link[4];
	/* For accept/reject */
#endif

#ifdef HAVE_U1
    complex u1tmp, loop, gftmp2; /* used by ../generic_u1/u1pot.c */
    /* the conjugate momentum for U(1) is a real(imaginary) number */
    Real mom_u1[4];
#endif
	/* The Kogut-Susskind phases, which have been absorbed into 
		the matrices.  Also the antiperiodic boundary conditions.  */
 	Real phase[4];

	/* 3 element complex vectors */
 	su3_vector phi[MAX_N_PSEUDO]; /* Gaussian random source, each pseudoferm */
 	su3_vector g_rand;	/* Gaussian random vector*/
	/* Use trick of combining xxx=D^adj D)^(-1) on even sites with
	   Dslash times this on odd sites when computing fermion force */
	
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
        /* store information about reunitarization */
        double RoS[4]; /* R^2/S^3 for cubic equation, normally =<1.0 */
        double gmin[4],gmax[4]; /* min/max eigenvalues of (V^+ V) */
        double denom[4]; /* denominator used for reunit force calculation */
        double unitW1[4]; /* deviation of W link from unitarity */
        int on_step_Y[4]; /* time step on which Y_phases were updated */
        int on_step_W[4]; /* time step on which W_norms were updated */
        int on_step_V[4]; /* time step on which V_dets were updated */
        Real phase_Y[4]; /* current Y matrix phase */
        Real phase_Y_previous[4]; /* previous Y matrix phase */
        su3_matrix Wlink[4];
        su3_matrix Wlink_previous[4];
        double Vdet[4]; /* abs of determinant of V matrix */
        double Xdet[4]; /* abs of determinant of X matrix before Naik */
        double XdetNaik[4]; /* abs of determinant of X matrix before Naik */
#endif /* HISQ_REUNITARIZATION_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */
} site;

/* End definition of site structure */

/* Definition of globals */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

/* The following are global scalars 
   beta is overall gauge coupling factor
   dyn_flavors are the number of flavors renormalizing the gauge action 
   u0 is tadpole improvement factor, perhaps (plaq/3)^(1/4)
*/
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  size_t volume;		/* volume of lattice = nx*ny*nz*nt */
#ifdef FIX_NODE_GEOM
EXTERN  int node_geometry[4];  /* Specifies fixed "nsquares" (i.e. 4D
			    hypercubes) for the compute nodes in each
			    coordinate direction.  Must be divisors of
			    the lattice dimensions */
#ifdef FIX_IONODE_GEOM
EXTERN int ionode_geometry[4]; /* Specifies fixed "nsquares" for I/O
			     partitions in each coordinate direction,
			     one I/O node for each square.  The I/O
			     node is at the origin of the square.
			     Must be divisors of the node_geometry. */
#endif
#endif
EXTERN  params param;
EXTERN	uint32_t iseed;		/* random number seed */
EXTERN  Real beta,u0;
EXTERN  int n_dyn_masses; // number of dynamical masses
EXTERN  Real dyn_mass[MAX_DYN_MASSES];
EXTERN  int dyn_flavors[MAX_DYN_MASSES];
#ifdef HAVE_U1
EXTERN  Real beta_u1;
/* Do we need a seperate u0_u1 for U1 links? */
EXTERN  int n_pseudo_charges; /* number of charges */
EXTERN  Real pseudo_charges[MAX_CHARGES];
EXTERN  int n_charges_uniq; /* number of unique charges */
EXTERN  Real charges_uniq[MAX_CHARGES];
EXTERN  Real current_charge_u1;
#endif
EXTERN	int warms,trajecs,steps,niter,nrestart,propinterval;
EXTERN  int niter_md[MAX_N_PSEUDO], niter_fa[MAX_N_PSEUDO], niter_gr[MAX_N_PSEUDO];
EXTERN  int prec_md[MAX_N_PSEUDO], prec_fa[MAX_N_PSEUDO], prec_gr[MAX_N_PSEUDO];
EXTERN  int prec_ff;
EXTERN	Real epsilon;
EXTERN	Real rsqmin_md[MAX_N_PSEUDO], rsqmin_fa[MAX_N_PSEUDO], rsqmin_gr[MAX_N_PSEUDO];
EXTERN  Real rsqprop;
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, RELOAD_BINARY,
			   RELOAD_CHECKPOINT, FRESH */
EXTERN	int saveflag;	/* do with lattice: FORGET, SAVE, SAVE_BINARY,
			   SAVE_CHECKPOINT */
/* Added by YL on 05/11/2016 */
#ifdef HAVE_U1
/* begin u(1) lattice:
 * fresh_u1, continue_u1 or reload_u1_ascii/serial/parallel
 */
EXTERN int start_u1flag;
/* end u(1) lattice:
 * save_u1_ascii/serial/parallel/scidac/parallel_scidac
 * or forget_u1
 */
EXTERN int save_u1flag;
EXTERN char start_u1file[MAXFILENAME];
EXTERN char save_u1file[MAXFILENAME];
#endif
/* end YL */

EXTERN	char rparamfile[MAXFILENAME];
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;

#ifdef HAVE_U1
EXTERN  Real g_splaq, g_tplaq;  /* global non-compact U1 plaquette measures */
EXTERN  Real *u1_A;             /* fundamental noncompact U1 field */
EXTERN  complex *link_u1;       /* noncompact U1 link exp(i * q * A) */
#ifdef HMC
EXTERN  Real *old_u1_A;
EXTERN  complex *old_link_u1;
/* For accept/reject */
#endif
EXTERN complex
*u1gf;           /* temporary global u1 field (complex vector potential) */
#endif

EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN	int total_iters;
EXTERN	int hisq_svd_counter;
EXTERN	int hisq_force_filter_counter;
EXTERN	int hyphisq_svd_counter;
EXTERN  int hypisq_force_filter_counter;
        /* source time, increment for it, and number of source slices */

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	size_t sites_on_node;		/* number of sites on this node */
EXTERN	size_t even_sites_on_node;	/* number of even sites on this node */
EXTERN	size_t odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

EXTERN gauge_file *startlat_p;
#ifdef HAVE_U1
EXTERN gauge_file *start_u1lat_p;
EXTERN gauge_file *save_u1lat_p;
#endif
EXTERN char hostname[128];

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 16
EXTERN char ** gen_pt[N_POINTERS];

/* Storage for definition of the quark action */
EXTERN fermion_links_t        *fn_links;
#ifdef HAVE_U1
/* EXTERN fermion_links_u1_t        *fn_links_u1; */
EXTERN imp_ferm_links_t *fn_combined;
#endif

EXTERN int n_pseudo;
EXTERN int max_rat_order;
EXTERN params_rhmc *rparam;
EXTERN int phases_in; /* 1 if KS and BC phases absorbed into matrices in site structure */

#ifdef MILC_GLOBAL_DEBUG
EXTERN int global_current_time_step;
#endif /* MILC_GLOBAL_DEBUG */
EXTERN int n_order_naik_total;
EXTERN int n_pseudo_naik[MAX_N_PSEUDO];
EXTERN int n_orders_naik[MAX_N_PSEUDO];
#ifdef HAVE_U1
EXTERN int n_pseudo_naik_charge[MAX_N_PSEUDO];
EXTERN int n_orders_naik_charge[MAX_N_PSEUDO];
EXTERN int n_pseudo_naik_charge_heavy[MAX_N_PSEUDO];
EXTERN int n_orders_naik_charge_heavy[MAX_N_PSEUDO];
#endif

/* For eigenpair calculation */
/* Not used for HMC */
EXTERN int Nvecs_tot;
EXTERN Real *eigVal; /* eigenvalues of D^dag D */
EXTERN su3_vector **eigVec; /* eigenvectors */

#endif /* _LATTICE_H */
