#ifndef _GENERIC_QUARK_TYPES_H
#define _GENERIC_QUARK_TYPES_H

#include "../include/macros.h"
#include "../include/precision.h"
#include "../include/int32type.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* Structures defining a generic quark source for both KS and Dirac fermions */

#define ALL_T_SLICES -1
#define MAXDESCRP 128
#define MAXSRCLABEL 8
#define MAXWEIGHTS 5

#ifdef HAVE_QIO
#include <qio.h>
#endif

/* For quark source and sink routines - both Wilson and KS */
/* The Weyl representation types are included for w_source_h */
enum source_type { 
  UNKNOWN = 0, 
  COMPLEX_FIELD_FILE, 
  COMPLEX_FIELD_FM_FILE, 
  COMPLEX_FIELD_STORE,
  CORNER_WALL, 
  CUTOFF_GAUSSIAN, 
  CUTOFF_GAUSSIAN_WEYL, 
  COVARIANT_GAUSSIAN,
  DERIV1,
  DERIV2_D,
  DERIV2_B,
  DERIV3_A,
  DIRAC_FIELD_FILE, 
  DIRAC_FIELD_FM_FILE, 
  DIRAC_FIELD_STORE,
  DIRAC_PROPAGATOR_FILE,
  EVEN_WALL, 
  EVENANDODD_WALL, 
  EVENMINUSODD_WALL, 
  FAT_COVARIANT_GAUSSIAN,
  FUNNYWALL1,
  FUNNYWALL2,
  GAUSSIAN, 
  HOPPING, 
  IDENTITY,
  POINT, 
  POINT_WEYL, 
  RANDOM_COLOR_WALL,
  RANDOM_CORNER_COLOR_WALL,
  ROTATE_3D,
  SPIN_TASTE,
  WAVEFUNCTION_FILE,
  VECTOR_FIELD_FILE, 
  VECTOR_FIELD_FM_FILE, 
  VECTOR_FIELD_STORE,
  VECTOR_PROPAGATOR_FILE
} ;

enum subset_type {
  FULL,
  HYPERCUBE
};

/* Header structure for a KS source in FNAL format */
typedef struct {
  int32type magic_number;
  int32type gmtime;
  int32type size_of_element;
  int32type elements_per_site;
  int32type dims[4];
  int32type site_order;
} ks_source_header;


/* File structure for a KS source in FNAL format */
typedef struct {
  ks_source_header *header;
  FILE *fp;
  char filename[MAXFILENAME];
  int byterevflag;
} ks_fm_source_file;

/* Structure defining a quark source operator */
struct qss_op_struct {
  int type;           /* operator type */
  char descrp[MAXDESCRP]; /* alpha description for most */
  char label[MAXSRCLABEL]; /* Abbreviation of description */
  Real a;             /* Lattice spacing for converting wave function file */
  Real d1;            /* Fermilab 3D rotation parameter */
  int dir1, dir2;     /* Directions for derivatives and hopping */
  int disp;           /* Stride for derivatives */
  Real weights[MAXWEIGHTS];  /* Weights for derivatives */
  Real eps_naik;      /* Naik epsilon for KS hopping operator */
  int dhop;           /* 0 for hop, 1 for 1st deriv of hop, 2 for 2nd */
  int iters;          /* iterations for covariant gaussian source */
  Real r0;            /* source size for gaussian, width for gauge invt  */
  int r_offset[4];    /* Coordinate offset for phases for some operators */
  int spin_taste;     /* For staggered fermions for some operators */
  char source_file[MAXFILENAME]; /* file name for some sources */
  struct qss_op_struct *op;   /* Next operation in the chain */
};

typedef struct qss_op_struct quark_source_sink_op;

/* Structure defining a staggered or Wilson (or clover) quark source */

typedef struct {
  int type;           /* source type */
  int orig_type;      /* original source type */
  int subset;         /* hypercube corners or full time slice */
  Real scale_fact;    /* scale factor */
  char descrp[MAXDESCRP];  /* alpha description for most */
  char label[MAXSRCLABEL]; /* Abbreviation of description */
  int ksource;        /* Counter for a list of wilson_vectors generated */
  int nsource;        /* Number of source wilson_vectors to be generated */
  int color;          /* Counter for source colors generated */
  int ncolor;         /* number of source su3_vectors to be generated*/
  int spin_snk;       /* Counter for KS propagators from a Dirac source */

  Real a;             /* Lattice spacing for converting wave function file */
  int x0,y0,z0,t0;    /* source coordinates for most */ 
  char source_file[MAXFILENAME]; /* file name for some sources */
  int sourceflag;      /* mode of reading or writing for some sources */
  char save_file[MAXFILENAME]; /* file name for saving the source */
  int savetype;        /* Type of file */
  int saveflag;           /* mode of writing */
  int source_file_initialized;
  int save_file_initialized;
  int mom[3];         /* insertion momentum for some sources */
  Real r0;            /* source size for gaussian, width for gauge invt  */
#ifdef HAVE_QIO
  QIO_Reader *infile;
  QIO_Writer *outfile;
#endif
  ks_fm_source_file *kssf;
  complex *c_src;      /* Pointer for complex source field storage */
  su3_vector *v_src;    /* su3_vector source for color walls */
  wilson_vector *wv_src; /* su3_vector source for color walls */
  quark_source_sink_op *op;   /* op need to create this 
				      source from parent */
  /* To be discontinued ... */
  int parity;         /* even or odd sites for w_source_h */
  int src_pointer ;   /* smearing function (for the moment, only
		         clover_finite_p_vary/create_wilson_source.c) */
  int wall_cutoff;    /* half size of box for w_source_h */

} quark_source;


/* Structure defining quark inversion parameters for most inverters */

typedef struct {
  int prec;           /* precision of the inversion 1 = single; 2 = double */
  int min;            /* minimum number of iterations (being phased out) */
  int max;            /* maximum number of iterations per restart */
  int nrestart;       /* maximum restarts */
  int parity;         /* EVEN, ODD, or EVENANDODD (for some inverters) */
  int start_flag;     /* 0: use a zero initial guess; 1: use dest */
  int nsrc;           /* Number of source vectors */
  Real resid;         /* desired residual - NOT SQUARED!
			 normalized as sqrt(r*r)/sqrt(src_e*src_e) */
  Real relresid;      /* desired relative residual - NOT SQUARED! */
  Real final_rsq;     /* Final true (absolute) residual. Norm'ed to (r*r)/(src*src) */
  Real final_relrsq;  /* Final relative residual. Same normalization. */
  Real size_r;        /* resulting cumulative residual. Same normalization. */
  Real size_relr;     /* resulting cumulative relative residual. Same normalization. */
  int converged;      /* returned 0 if not converged; 1 if converged */
  int  final_iters;
  int  final_restart;
                      /* Add further parameters as needed...  */
} quark_invert_control;

void report_status(quark_invert_control *qic);

/* Structures required for specific inverters */

/* Structure defining parameters of Dirac matrix for clover inversion */
/* To be passed through to inverter. */
typedef struct {
  Real Kappa;        /* hopping */
  Real Clov_c;       /* Perturbative clover coeff */
  Real U0;           /* Tadpole correction to Clov_c */
} dirac_clover_param;

/* Same for Wilson case */
typedef struct {
  Real Kappa;        /* hopping */
} dirac_wilson_param;

/* Same for plain KS case */
typedef struct {
  Real mass;
  Real offset;
  int naik_term_epsilon_index;
  Real naik_term_epsilon;

} ks_param;

#endif /* _GENERIC_QUARK_TYPES_H */


