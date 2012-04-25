#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice_form.h ********************************/
/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice. */

/* Modifications
   C. McNeile 1997 Original version
   C. DeTar 5/24/97 Added quark_zonked_rot to site structure 
   */

#include "defines.h"
#include "../include/generic_wilson.h" 
#include "../include/generic_quark_types.h"
#include "../include/random.h"
#include "../include/macros.h"
#include "../include/io_lat.h"    /* For gauge_file */

/* Begin definition of site structure */

#include "../include/su3.h"
#include "../include/complex.h"

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


	/* wilson half vector (temporary used in dslash) */
	/*LEAN*half_wilson_vector htmp[8];**/
	half_wilson_vector htmp[2];

     /****** variables for the hqet--> light  code ***/

        /* quark propagator inversion workspace ****/
	wilson_vector psi ;
	wilson_vector mp ;

	spin_wilson_vector quark_zonked ;
	spin_wilson_vector quark_zonked_rot ;
	spin_wilson_vector quark_sequential ;
	spin_wilson_vector quark_spectate ;

#ifdef BICG_CLOVER
	wilson_vector sss;	/* internal biconjugate gradient vector */
	wilson_vector tmpb;	/* auxiliary for other internal biCG vectors */
	wilson_vector tmp ;	/* auxiliary for other internal biCG vectors */
#endif


	/** work space for the hqet inversion ***/
	su3_matrix tempvec[4];


        /* hqet propagator ******/
	su3_matrix heavy_prop ; 

	/* The smearing functions for the sequential source  ****/
#define MAXMOM   20
	complex seq_smear_func[ MAXMOM ] ;
	complex seq_smear_func_fft[ MAXMOM ] ;

	/* FFT work space *******/
	complex seq_smear_fft1[ MAXMOM ] ;
	complex seq_smear_fft2[ MAXMOM ] ;

	/** Work space for the smearing of the sequential source  ******/
        wilson_vector fft_one;
        wilson_vector fft_two;

        /** The source for the propagator inversions ***/
#define chi fft_one
#define src_store  fft_two
} site;


/* End definition of site structure */

/* Definition of globals */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

#ifdef CONTROL
int n_spins = 4;                /* Number of spins generated */
int source_loc[4] = {0,0,0,0};  /* Source location for all propagators */
int spins[4] = {0,1,2,3};       /* List of spins generated */
#else
extern int n_spins;
extern int source_loc[4];
extern int spins[4];
#endif
/***  #define EPS 1.0e-6   MORE WORK *****/

enum quark_type_choices { STANDARD_QUARK , STATIC_LIGHT_QUARK  } ; 
#define UNITMATRIX -100


/***** flag that controls whether the source gamma matrix
       representation is flipped ****/
/**#define FLIP_Q_SOURCE_REP  1  ***/

/* The following are global scalars */
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  int volume;	        /* volume of lattice = nx*ny*nz*nt */

EXTERN	Real clov_c,u0;
EXTERN  Real byterevReal;  /* flag to bytereverse lattice */

EXTERN	char startfile[MAXFILENAME] ; /*** file containing the gauge configurations ***/
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN  int fixflag;  /* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */

EXTERN  int no_spectator ;  /** number of spectator kappa values ***/
EXTERN  int no_zonked_light ; /*** the number of light kappa values of
				zonked quarks ***/

/**** The maximum number of kappa values ***/
#define MAX_KAPPA 10
EXTERN Real kappa_zonked_light[MAX_KAPPA] ; /** Kappa values for the
					       light zonked quark **/
EXTERN Real kappa_spectator[MAX_KAPPA] ; /** Kappa values for the
					    spectator quarks inversion
					    **/

EXTERN int startflag_spectator[MAX_KAPPA] ;  /** type of IO read for
				      the spectator quark **/
EXTERN int startflag_zonked[MAX_KAPPA] ;  /** type of IO read for the
					      zonked quark **/
EXTERN Real resid_spectator;
EXTERN Real resid_zonked,r0_zonked[MAX_KAPPA];

EXTERN	char qfile_spectator[MAX_KAPPA][MAXFILENAME];  
EXTERN	char qfile_zonked[MAX_KAPPA][MAXFILENAME];

/** The maximum number of velocities ****/
#define MAXVEL  8  
EXTERN	char hqet_smear_file[MAXVEL ][MAXFILENAME]; /** File containing the
					     sources for the
					     sequential inversion
					     ****/


EXTERN	char heavy_light_out[MAXFILENAME] ; /** File to write the heavy -->
				     light form factors to **/
EXTERN	char twopt_out[MAXFILENAME] ;       /** File to write the two point
				     functions ***/
EXTERN	char seq_out[MAXFILENAME] ;       /** File to write the sequential two
				   point functions ***/

EXTERN int no_q_values ;
EXTERN int quark_type  ;

EXTERN int verbose_flag ; /*** flag controllling the amount of debug
			    information to print ***/

/*** THe number of momentum values used in the code ****/
EXTERN int q_momstore[MAXMOM][3] ;

EXTERN  Real velocity[MAXVEL][4] ;  /*** The store of velocities ******/
EXTERN  int novel ;   /** The number of velocities   ****/

EXTERN int tf ; /** Fixed timeslice in the sequential source ***/

EXTERN	Real width,source_r0 ;  /*** mORE work ***/


/*** stuff for the inversion check on the quark propagators *****/
EXTERN  int source_parity;
EXTERN  int wall_cutoff,wall_separation;
EXTERN  Real rsqmin,rsqprop ;
EXTERN  int niter_spectator,niter_zonked,
  nrestart_spectator,nrestart_zonked ; 
EXTERN  Real kappa; /*nickname for various kappa values in the loops ***/
EXTERN  int wallflag;

EXTERN int layout_flag ;
enum layout_flag_set { HYPER_GRAY_EVENFIRST = 100 } ;

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN  int	sites_on_node;		/* number of sites on this node */
EXTERN  int	even_sites_on_node;	/* number of even sites on this node */
EXTERN  int	odd_sites_on_node;	/* number of odd sites on this node */
EXTERN  int	number_of_nodes;	/* number of nodes in use */
EXTERN  int  this_node;		/* node number of this node */
EXTERN Real *wall_template; /* pointer to array of wall heights */

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;

EXTERN quark_source wqs_spectator[MAX_KAPPA];
EXTERN quark_source wqs_zonked_light[MAX_KAPPA];
EXTERN quark_source wqstmp;  /* Temporary */

EXTERN quark_invert_control qic_zonked_light;
EXTERN quark_invert_control qic_spectator;

EXTERN dirac_clover_param dcp;
EXTERN dirac_wilson_param dwp;

EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8	/* Number of generic pointers */
EXTERN char ** gen_pt[N_POINTERS];


/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

#endif /* _LATTICE_H */
