#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD version 7
   clover_hybrids program
   This file defines global scalars and the fields in the lattice. */

#include "defines.h" /* For SITERAND */
#include "../include/generic_quark_types.h"
#include "../include/macros.h"
#include "../include/random.h"
#include "../include/generic_clover.h" /* For clover */

/* Begin definition of site structure */

#include "../include/su3.h"
#include "../include/random.h"   /* For double_prn */

/* The lattice is an array of sites.  */
typedef struct  {
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

	/* Use unions so smeared gauge fields overlap wilson vectors.
	   Programmer must use only one at a time!!  In particular, this
	   means that you must do all your smearing, and compute the
	   field strength, before doing any conjugate gradients, since
	   the CG will overwrite the smeared links. */
	union {
	    wilson_vector Uvec[5];	/* cg vectors */
	    su3_matrix Umat[6];	/* smearlinks, hyb_tempmat1 and hyb_tempmat2 */
	} U;
	union {
	    wilson_vector Vvec[3];	/* */
	    su3_matrix Vmat[4];	/* */
	} V;

	/* Wilson complex vectors */
#define	      G_RAND  U.Uvec[0]	/* gaussian random vector */
#define	      PSI  U.Uvec[1]	/* solution vector */
#define	      CHI  U.Uvec[2]	/* source vector */
#define	      MP   U.Uvec[3]	/* another CG vector */
	wilson_vector ttt;	/* internal conjugate gradient vector */
#define		      vtmp  U.Uvec[4]	/* more cg vectors for bicongrad */
	wilson_vector sss;	/* more cg vectors for bicongrad */

        half_wilson_vector htmp[2];

#ifdef SMEAR
	/* smeared links, and temporary matrix */
#define		   smearlink  U.Umat	/* su3_matrix[4] */
#define		   templink   V.Vmat	/* su3_matrix[4] */
#define		   hyb_tempmat1 U.Umat[4]	/* su3_matrix */
#define		   hyb_tempmat2 U.Umat[5]	/* su3_matrix */
#endif
        /* components of F_mu_nu, use defines to give  components names */
        su3_matrix field_strength[6];
        /* quark and antiquark propagators */
#define		quark_source V.Vvec[0]	/* wilson_vector */
#define		quark_prop V.Vvec[1]	/* wilson_vector */
#define		anti_prop V.Vvec[2]	/* wilson_vector */
  /* To allow cg restarts.  Saves source... CD */
#define         quark_save U.Uvec[0]    /* wilson_vector */


/**** clover stuff ***/
        wilson_vector tmp;      /* another temporary CG vector */
         wilson_vector tmpb;    /* auxiliary for other internal biCG vectors */
/**** end of clover stuff ***/

	/** storage space for the baryon correlators ***/
	wilson_propagator quark_store ;

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
EXTERN	int niter;
EXTERN	int source_start, source_inc, n_sources;
	/* source time, increment for it, and number of source slices */
EXTERN	Real rsqprop,beta,kappa;
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN	int startflag;	/* beginning lattice flag */
EXTERN  int saveflag;	/* end lattice flag */
EXTERN	int total_iters;
EXTERN  int wot_src ;     /*** source for the quark propagator inverter ***/

EXTERN	Real clov_c,u0;
EXTERN  int fixflag;  /* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */
EXTERN  int boundary_flag  ;
EXTERN int verbose_flag ; /*** flag controlling the amount of debug information to print ***/

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

/*** factors controlling the smearing level ***/
EXTERN  Real  space_simple_weight  ;
EXTERN  Real  space_norm_factor   ;
EXTERN  Real  time_simple_weight ;
EXTERN  Real  time_norm_factor   ;
EXTERN  int  smearing_level ;


/*** list of which hybrid operators to calculate ***/

EXTERN int oper_PION_SOURCE  ;
EXTERN int oper_PION2_SOURCE  ;
EXTERN int oper_RHO_SOURCE   ; 
EXTERN int oper_RHO2_SOURCE ;
EXTERN int oper_A1P_SOURCE   ;
EXTERN int oper_A1_SOURCE   ;
EXTERN int oper_ZEROMP_SOURCE ;
EXTERN int oper_ZEROPM_SOURCE ;
EXTERN int oper_ZEROPMP_SOURCE ;
EXTERN int oper_ZEROPMB_SOURCE ;
EXTERN int oper_ZEROMM_SOURCE ; 
EXTERN int oper_ZEROMMP_SOURCE ;
EXTERN int oper_ONEMP_SOURCE  ;
EXTERN int oper_ONEMP2_SOURCE  ;
EXTERN int oper_ONEMM_SOURCE  ;
EXTERN int oper_ONEPP_SOURCE  ;
EXTERN int oper_QQQQ_SOURCE   ;

enum operator_choices { CALCULATE = 100 , DO_NOT_CALCULATE = 101 } ; 

/** The minimum number of iterations should not be used in this code **/
#ifdef CONTROL
int min_iters = -1 ; 
#else
extern  int min_iters ; 
#endif

EXTERN quark_invert_control qic;
EXTERN dirac_clover_param dcp;

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

/* defines for index on field_strength */
#define FS_XY 0
#define FS_XZ 1
#define FS_YZ 2
#define FS_XT 3
#define FS_YT 4
#define FS_ZT 5

 
/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8	/* Number of generic pointers */
/* NEED 8 WHEN GAUGEFIXING */
EXTERN char ** gen_pt[N_POINTERS];

/* Storage for the clover term */
EXTERN clover *gen_clov;

#endif /* _LATTICE_H */
