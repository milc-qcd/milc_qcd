#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/
/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/generic_wilson.h" 
#include "../include/generic_quark_types.h"
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/random.h"    /* For double_prn */
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
#define MAXHTMP 8
	/*LEAN*half_wilson_vector htmp[8];**/
	half_wilson_vector htmp[MAXHTMP];

     /****** variables for the propagating form factor code ***/

        /* the quark propagator  ****/
	wilson_vector mp ;

        spin_wilson_vector quark_zonked_heavy[MAX_ZONKED_HEAVY];
        spin_wilson_vector quark_zonked_light[MAX_ZONKED_LIGHT];
	spin_wilson_vector quark_sequential ;
	spin_wilson_vector quark_spectator ;
	spin_wilson_vector quark_rot ; 

#ifdef LARGE_MEMORY_SMEAR
#define large_fft_one quark_rot
	spin_wilson_vector large_fft_two;
#endif

#ifdef CLOVER
	wilson_vector sss;	/* internal biconjugate gradient vector */
	wilson_vector tmpb;	/* auxiliary for other internal biCG vectors */
	wilson_vector tmp ;	/* auxiliary for other internal biCG vectors */

#endif

	/* The smearing functions for the sequential source  ****/
	complex seq_smear_func[ MAXPMOM ] ;

	/* smearing functions in coord space for heavy quarks ****/
        /** NOTE: We use only the [0] element for now **/
	complex heavy_smear_func[MAX_ZONKED_HEAVY] ; 

	/* smearing functions in mom space for heavy quarks ****/
        /** NOTE: We use only the [0] element for now **/
	complex heavy_smear_func_mom[MAX_ZONKED_HEAVY] ; 

	/* smearing functions in coord space for light quarks ****/
        /** NOTE: Unused for now **/
	complex light_smear_func[MAX_ZONKED_LIGHT] ; 

	/* smearing functions in mom space for light quarks ****/
        /** NOTE: Unused for now **/
	complex light_smear_func_mom[MAX_ZONKED_LIGHT] ; 

	/** Work space for the Wilson line load ******/
        /* CAUTION!! There is an implicit union going on
	   here - This needs to be fixed. -see QTMP below */
        wilson_vector psi;
        wilson_vector fft_one;
        wilson_vector fft_two;
        su3_matrix staple;   /* so we can calculate the plaquette */

        /***  used for quark inversion ***/
#define chi fft_one
#define src_store  fft_two

        /***  used in meson_cont_mom_lean2  ***/

#define QTMP psi     /* Space for FT phases used in meson_cont_mom_lean2 */
#define DIMQTMP 45   /* Number of complex values that fit in QTMP        */

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
EXTERN int wallflag;

/* The following are global scalars */
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  int volume;	/* volume of lattice = nx*ny*nz*nt */

EXTERN int verbose_flag ; /*** flag controllling the amount of debug
			    information to print ***/
EXTERN Real kappa ; 

EXTERN	Real clov_c,u0;
EXTERN  Real byterevReal;  /* flag to bytereverse lattice */

EXTERN	char startfile[MAXFILENAME] ; /*** file containing the gauge configurations ***/
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN  int fixflag;

EXTERN  int no_spectator ;  /** number of spectator kappa values ***/
EXTERN  int no_zonked_light ; /*** the number of light kappa values of zonked quarks ***/
EXTERN  int no_zonked_heavy ; /*** the number of light kappa values of zonked quarks ***/
EXTERN  int no_sequential ; /** number of sequential kappa values to use ***/

EXTERN Real kappa_spectator[MAX_KAPPA] ; /** Kappa values for the spectator quarks inversion  **/
EXTERN Real kappa_zonked_light[MAX_KAPPA] ; /** Kappa values for the light zonked quark **/
EXTERN Real kappa_zonked_heavy[MAX_KAPPA] ; /** Kappa values for the heavy zonked quark **/
EXTERN Real kappa_sequential[MAX_KAPPA] ; /** Kappa values for the sequential inversion  **/

EXTERN int startflag_spectator[MAX_KAPPA]   ;  /** type of IO read for the spectator quark ***/
EXTERN int startflag_zonked_light[MAX_KAPPA]   ;  /** type of IO read for the zonked quark ***/
EXTERN int startflag_zonked_heavy[MAX_KAPPA]   ;  /** type of IO read for the zonked quark ***/
EXTERN int saveflag_zonked_heavy[MAX_KAPPA]   ;  /** type of IO read for the zonked quark ***/
EXTERN int saveflag_zonked_light_ssink   ;  /** type of IO read for the light zonked quark ***/
EXTERN int reloadflag_zonked_light_ssink   ;  /** type of IO read for the light zonked quark ***/
EXTERN int saveflag_zonked_heavy_ssink   ;  /** type of IO read for the heavy zonked quark ***/
EXTERN int reloadflag_zonked_heavy_ssink  ;  /** type of IO read for the heavy zonked quark ***/

						
EXTERN Real resid_spectator;
EXTERN Real resid_zonked_light;
EXTERN Real resid_zonked_heavy;

EXTERN	char qfile_spectator[MAX_KAPPA][MAXFILENAME];  
EXTERN	char qfile_zonked_light[MAX_KAPPA][MAXFILENAME];
EXTERN	char qfile_zonked_light_ssink[MAX_KAPPA][MAXFILENAME]; /** generated inside the code **/
EXTERN	char qfile_suffix_zonked_light[MAXFILENAME];
EXTERN	char qfile_zonked_heavy[MAX_KAPPA][MAXFILENAME]; /** generated inside the code **/
EXTERN	char qfile_zonked_heavy_ssink[MAX_KAPPA][MAXFILENAME]; /** generated inside the code **/
EXTERN	char qfile_suffix_zonked_heavy[MAXFILENAME];

EXTERN	char seq_smear_file[MAXPMOM][MAXFILENAME]; /** File containing the  sources for the sequential inversion ****/


EXTERN	char filename_HH3[MAXFILENAME] ; /** File to write the heavy --> heavy form factors to ***/
EXTERN	char filename_HL3[MAXFILENAME] ; /** File to write the heavy --> light form factors to ***/

EXTERN int saveflag_HH3;       /* Flags to specify output type */
EXTERN int saveflag_HL3;       /* choices are FORGET, SAVE_ASCII */

EXTERN	char filename_HH2_GL[MAXFILENAME] ; /** File to write the heavy-heavy two point functions  to ***/
EXTERN	char filename_LL2_GG[MAXFILENAME] ; /** File to write the light_light two point functions  to ***/
EXTERN	char filename_HL2_GG[MAXFILENAME] ; /** File to write the heavy-light two point functions  to ***/
EXTERN	char filename_HL2_GE[MAXFILENAME] ; /** File to write the heavy-light two point functions  to ***/
EXTERN  char filename_HL2_GL[MAXFILENAME] ; /** File to write the heavy-light two point functions to ***/
  				  
EXTERN int saveflag_HH2_GL;       /* correlator file type */
EXTERN int saveflag_LL2_GG;
EXTERN int saveflag_HL2_GG;
EXTERN int saveflag_HL2_GE;
EXTERN int saveflag_HL2_GL;

EXTERN int no_p_values ;
EXTERN int no_q_values ;
EXTERN int no_k_values ;

EXTERN int p_momstore[MAXPMOM][3] ;
EXTERN int q_momstore[MAXMOM][3] ;
EXTERN int k_momstore[MAXMOM][3] ;

EXTERN int tf ; /** Fixed timeslice in the sequentiaol source ***/

EXTERN  char walldescrp[3][30];  /* Describes type of source */

EXTERN  int source_x,source_y,source_z,source_t; /*source points */

EXTERN	int niter_spectator,niter_zonked_light,niter_zonked_heavy,
  nrestart_spectator,nrestart_zonked_light,nrestart_zonked_heavy ;


/*** io look up table for saving the temporary ***/
/**typedef struct { int spin ; int color ; } lookup_w_io   ;
   EXTERN lookup_w_io first_loc[4] , second_loc[4] ; **//** defined in the code **/

EXTERN int *w_meson_store_t;  /* Map time slice to index */
EXTERN int *w_meson_my_t;     /* Inverse map - index to time slice */
EXTERN int w_meson_nstore;    /* Max index stored on this node */
EXTERN complex *w_meson_corr; /* Storage for correlator        */

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN  int	sites_on_node;		/* number of sites on this node */
EXTERN  int	even_sites_on_node;	/* number of even sites on this node */
EXTERN  int	odd_sites_on_node;	/* number of odd sites on this node */
EXTERN  int	number_of_nodes;	/* number of nodes in use */
EXTERN  int  this_node;		/* node number of this node */

EXTERN quark_source wqs_spectator[MAX_KAPPA];
EXTERN quark_source wqs_zonked_light[MAX_KAPPA];
EXTERN quark_source wqs_zonked_light_tmp[MAX_KAPPA];
EXTERN quark_source wqs_zonked_heavy[MAX_KAPPA];
EXTERN quark_source wqs_zonked_heavy_tmp[MAX_KAPPA];
EXTERN quark_source wqstmp;  /* Temporary */

EXTERN quark_invert_control qic_zonked_heavy;
EXTERN quark_invert_control qic_zonked_light;
EXTERN quark_invert_control qic_sequential;
EXTERN quark_invert_control qic_spectator;

EXTERN int inverter_type_zonked_heavy[MAX_KAPPA];
EXTERN int inverter_type_sequential[MAX_KAPPA];

EXTERN dirac_clover_param dcp;
EXTERN dirac_wilson_param dwp;

EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8	/* Number of generic pointers */
EXTERN char ** gen_pt[N_POINTERS];

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

#define IF_MASTER  if(this_node==0)

#endif /* _LATTICE_H */
