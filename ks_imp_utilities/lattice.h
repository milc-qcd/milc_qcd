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
     Includes storage for "fat links"  (fatlink[4])
*/

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/random.h"    /* For double_prn */
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */

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

	/* antihermitian momentum matrices in each direction */
 	anti_hermitmat mom[4];
#ifdef FERMION_FORCE
        su3_matrix ansmom[4];
#endif

	/* The Kogut-Susskind phases, which have been absorbed into 
		the matrices.  Also the antiperiodic boundary conditions.  */
 	Real phase[4];

	/* 3 element complex vectors */
 	su3_vector phi;	        /* Gaussian random source vector */
 	su3_vector resid;	/* conjugate gradient residual vector */
 	su3_vector cg_p;	/* conjugate gradient change vector */
 	su3_vector xxx;	        /* solution vector = Kinverse * phi */
 	su3_vector ttt;		/* temporary vector, for K*ppp */
 	su3_vector g_rand;	/* Gaussian random vector*/
	/* Use trick of combining xxx=D^adj D)^(-1) on even sites with
	   Dslash times this on odd sites when computing fermion force */

	/* temporary vectors and matrices */
	su3_vector tempvec[4];	/* One for each direction */
#ifdef FN
	su3_vector templongvec[4];	/* One for each direction */
        su3_vector templongv1;
#endif
	su3_matrix tempmat1,staple;
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
   mass1 and mass2 are quark masses
   u0 is tadpole improvement factor, perhaps (plaq/3)^(1/4)
*/
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  int volume;		/* volume of lattice = nx*ny*nz*nt */
EXTERN	int iseed;		/* random number seed */
EXTERN	int niter;
EXTERN  Real beta;
EXTERN  Real mass,u0;
EXTERN	Real rsqprop;
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, RELOAD_BINARY,
			   RELOAD_CHECKPOINT, FRESH */
EXTERN	int saveflag;	/* do with lattice: FORGET, SAVE, SAVE_BINARY,
			   SAVE_CHECKPOINT */
EXTERN  int savelongflag, savefatflag;
EXTERN  int srcflag, ansflag;
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN  char savelongfile[MAXFILENAME],savefatfile[MAXFILENAME];
EXTERN  char srcfile[MAXFILENAME],ansfile[MAXFILENAME];
EXTERN  int inverttype;
EXTERN	int total_iters;
EXTERN  int phases_in; /* 1 if KS and BC phases absorbed into matrices */

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

#ifdef HAVE_QDP
#include <qdp.h>

EXTERN QDP_ColorMatrix **fatlinks, **longlinks, *implinks[8];
EXTERN QDP_Shift neighbor3[4];
EXTERN QDP_Shift shiftdirs[8];
EXTERN QDP_ShiftDir shiftfwd[8], shiftbck[8];

#endif

#endif /* _LATTICE_H */
