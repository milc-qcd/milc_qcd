#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD version 6
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
#include "../include/random.h"    /* For double_prn */
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */
#ifdef HAVE_QDP
#include <qdp.h>
#endif

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

/* ------------------------------------------------------------ */
/*   Now come the physical fields, program dependent            */
/* ------------------------------------------------------------ */
	/* gauge field */
	su3_matrix link[4];	/* the fundamental field */
        su3_matrix rgt;  /* gauge transformation */
        su3_matrix sm_link[4];  /* the smeared field */
        su3_matrix stapleg;
        su3_vector tempvecg;       
#if  defined(HYBRIDS)
	su3_matrix longlink[4];	/* three link straight paths */
#endif

#ifdef HMC_ALGORITHM
 	su3_matrix old_link[4];
	/* For accept/reject */
#endif

	/* antihermitian momentum matrices in each direction */
 	anti_hermitmat mom[4];

	/* The Kogut-Susskind phases, which have been absorbed into 
		the matrices.  Also the antiperiodic boundary conditions.  */
 	Real phase[4];

	/* 3 element complex vectors */
 	su3_vector phi;		/* Gaussian random source vector */
 	su3_vector xxx;		/* solution vector = Kinverse * phi */
#ifdef PHI_ALGORITHM
 	su3_vector old_xxx;	/* For predicting next xxx */
#endif
 	su3_vector resid;	/* conjugate gradient residual vector */
 	su3_vector phi2;	/* Gaussian random source vector, mass2 */
 	su3_vector cg_p;	/* conjugate gradient change vector */
 	su3_vector ttt;		/* temporary vector, for K*ppp */
 	su3_vector g_rand;	/* Gaussian random vector*/
	/* Use trick of combining xxx=D^adj D)^(-1) on even sites with
	   Dslash times this on odd sites when computing fermion force */
	
#ifdef HYBRIDS
        su3_matrix field_strength[6];
#endif
	/* temporary vectors and matrices */
	su3_vector tempvec[4];	/* One for each direction */
        su3_matrix s_link,s_link_f,t_link_f,diag,test[4];
#ifdef FN
	su3_vector templongvec[4];	/* One for each direction */
        su3_vector templongv1;
#endif
#ifdef DM_DU0
	su3_vector dMdu_x;	/* temp vector for <psi-bar(dM/du0)psi>
				   calculation in 'f_meas.c' */
#endif
#ifdef NPBP_REPS
 	su3_vector M_inv;	/* temp vector for M^{-1} g_rand */
#endif
#ifdef CHEM_POT
 	su3_vector dM_M_inv;	/* temp vector for dM/dmu M^{-1} g_rand */
#endif
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
   mass, mass1 and mass2 are quark masses
   u0 is tadpole improvement factor, perhaps (plaq/3)^(1/4)
*/
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  int volume;		/* volume of lattice = nx*ny*nz*nt */
EXTERN	int iseed;		/* random number seed */
EXTERN	int warms,trajecs,steps,niter,propinterval;
EXTERN  int nflavors;
EXTERN	Real epsilon;
EXTERN  Real beta,u0;
EXTERN  Real mass;
EXTERN	Real rsqmin,rsqprop;
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, RELOAD_BINARY,
			   RELOAD_CHECKPOINT, FRESH */
EXTERN	int saveflag;	/* do with lattice: FORGET, SAVE, SAVE_BINARY,
			   SAVE_CHECKPOINT */
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN	int total_iters;
EXTERN  int phases_in; /* 1 if KS and BC phases absorbed into matrices */
EXTERN  int source_start, source_inc, n_sources;
        /* source time, increment for it, and number of source slices */

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

/* Flags to help us remember what we are doing */
EXTERN int valid_longlinks;
EXTERN int valid_fatlinks;

EXTERN gauge_file *startlat_p;

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;


/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN site *lattice;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 16
EXTERN char ** gen_pt[N_POINTERS];

/* field major storage DON't FORGET to MALLOC somewhere */
EXTERN su3_matrix *t_longlink;
EXTERN su3_matrix *t_fatlink;
#ifdef DBLSTORE_FN
EXTERN su3_matrix *t_longbacklink;
EXTERN su3_matrix *t_fatbacklink;
#endif
#ifdef DM_DU0
EXTERN su3_matrix *t_dfatlink_du0;
#endif

#ifdef HAVE_QDP
EXTERN QDP_ColorMatrix **fatlinks, **longlinks, *implinks[8];
EXTERN QDP_Shift neighbor3[4];
EXTERN QDP_Shift shiftdirs[8];
EXTERN QDP_ShiftDir shiftfwd[8], shiftbck[8];
#endif
#define smear_fac 2.5
#endif /* _LATTICE_H */
