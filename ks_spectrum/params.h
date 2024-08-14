#ifndef _PARAMS_H
#define _PARAMS_H

#include "defines.h"
#include "../include/macros.h"  /* For MAXFILENAME */
#include "../include/gb_ops.h"
#include "../include/generic_quark_types.h"
#include "../include/generic_ks.h" /* For quark_source */
#include "../include/generic_wilson.h"  /* For quark_source */
#include "../include/gammatypes.h"
#include "../include/imp_ferm_links.h"

#define MAX_MASS_PBP 8
#define MAX_SOURCE 128
#define MAX_SET 256
#define MAX_PROP 64
#define MAX_QK 512
#define MAX_COMBO 8
#define MAX_PAIR 5000
#define MAX_TRIPLET 64
#define MAX_QKPAIR_LABEL 64
#define MAX_MESON 32
#define MAX_SPECTRUM_REQUEST 512
#define MAX_MESON_LABEL 64
#define MAX_BARYON_LABEL 64
#define MAX_MESON_MOMENTUM 100
#define MAX_MOM_LABEL 32
#define MAX_CORR 256
#define STATIC_TYPE 0
#define KS_TYPE 1
#define PROP_TYPE 0
#define QUARK_TYPE 1
#define COMBO_TYPE 2
#define BASE_SOURCE_PARENT -1

#ifdef GB_BARYON
/* specific to golterman-bailey baryon 2-pts */
#define MAX_OCTET 64
#endif

enum check_type { CHECK_NO,  CHECK_YES, CHECK_SOURCE_ONLY };

enum set_type { MULTIMASS_SET, MULTISOURCE_SET, MULTICOLORSOURCE_SET, SINGLES_SET };

/* structure for passing simulation parameters to each node */
typedef struct {
  int stopflag;   /* 1 if it is time to stop */
  /* INITIALIZATION PARAMETERS */
  int nx,ny,nz,nt;	/* lattice dimensions */
#ifdef FIX_NODE_GEOM
  int node_geometry[4];  /* Specifies fixed "nsquares" (i.e. 4D
			    hypercubes) for the nodes in each
			    coordinate direction.  Must be divisors of
			    the lattice dimension */
#ifdef FIX_SUBNODE_GEOM
  int subnode_geometry[4];  /* Specifies fixed "nsubsquares" (i.e. 4D
			    hypercubes) for the PE ranks on each node in each
			    coordinate direction.  Must be divisors of
			    the node sublattice dimensions -- that is
			    full lattice dims divided by node_geometry */
#endif
#ifdef FIX_IONODE_GEOM
  int ionode_geometry[4]; /* Specifies fixed "nsquares" for I/O
			     partitions in each coordinate direction,
			     one I/O node for each square.  The I/O
			     node is at the origin of the square.
			     Must be divisors of the node_geometry. */
#endif
#endif
  uint32_t iseed;
  char job_id[MAXFILENAME]; /* Usually encoded by scripts */

  /*  REPEATING BLOCK */
  int startflag;	/* what to do for beginning lattice */
  int start_u1flag;	/* what to do for beginning u(1) lattice */
  int start_ape_flag;   /* what to do for input APE links */
  Real u0;
  int coord_origin[4];  /* Origin of coordinates for KS phases and time_bc */
  int time_bc;          /* 0 for antiperiodic, 1 for periodic */
  int fixflag;          /* whether to gauge fix */
  int saveflag;	        /* what to do for saving lattice */
  int save_u1flag;	/* what to do with ending u(1) lattice */
  int save_ape_flag;	/* how to save the APE links */
  Real staple_weight;
  int ape_iter;
#if EIGMODE == EIGCG
  eigcg_params eigcgp; /* parameters for eigCG */
#endif
  int n_naiks;                /* For HISQ links */
  double eps_naik[MAX_NAIK];  /* For HISQ links */
  char ks_eigen_startfile[MAXFILENAME]; /* KS eigenvector file to be loaded */
  char ks_eigen_savefile[MAXFILENAME]; /* KS eigenvector file to be saved */
  int ks_eigen_startflag; /* what to do for beginning eigenvectors */
  int ks_eigen_saveflag; /* what to do for ending eigenvectors */
  ks_eigen_param eigen_param; /* Parameters for eigensolver */
  int num_pbp_masses;   /* Number of masses for pbp calculation */
  quark_invert_control qic_pbp[MAX_MASS_PBP];
  int prec_pbp;         /* Precision of the pbp calculation (1 or 2) */
  int npbp_reps;     /* Number of random sources for pbp calculation */
  ks_param ksp_pbp[MAX_MASS_PBP];
  Real charge_pbp[MAX_MASS_PBP];
  /* Sources */
  int num_base_source;  /* Number of base sources */
  quark_source src_qs[MAX_SOURCE]; /* source parameters */
  //  quark_source base_src_qs[MAX_SOURCE];
  int num_modified_source;       /* Number of modified sources */
  quark_source_sink_op src_qs_op[MAX_SOURCE];
  int parent_source[MAX_SOURCE];      /* base_source or source index */
  /* Multimass or multisource sets */
  int num_set;  /* number of sets */
  enum set_type set_type[MAX_SET];    /* multimass or multisource */
  enum inv_type inv_type[MAX_SET];    /* inverter type MG or CG */
  Real charge[MAX_SET];     /* charge for propagators in the set */
  char charge_label[MAX_SET][32];  /* for correlator label */
  int num_prop[MAX_SET]; /* number of propagators in a set */
  int begin_prop[MAX_SET]; /* index of beginning propagator in this set */
  int end_prop[MAX_SET]; /* index of ending propagator in this set */
  /* Propagators */
  int prop_type[MAX_PROP]; /* 0 static 1 KS */
  int set[MAX_PROP];  /* The set to which the propagator belongs */
  int startflag_ks[MAX_PROP];	/* what to do for beginning KS prop */
  int saveflag_ks[MAX_PROP];	/* what to do for saving KS prop */
  int source[MAX_PROP];      /* index of source for this prop */
  char mass_label[MAX_PROP][32]; /* mass label for this prop */
  ks_param ksp[MAX_PROP];         /* propagator parameters for this prop */
  quark_invert_control qic[MAX_PROP];
  enum check_type check[MAX_PROP];         /* True -> run the inverter */
  Real bdry_phase[MAX_PROP][4];      /* For twisted boundary conditions */
  char startfile_ks[MAX_PROP][MAXFILENAME];
  char savefile_ks[MAX_PROP][MAXFILENAME];
  /* Quarks */
  int num_qk;	                     /* number of quarks */
  quark_source_sink_op snk_qs_op[MAX_QK];
  int parent_type[MAX_QK];           /* propagator type: quark or propagator */
  int prop_for_qk[MAX_QK];           /* Propagator or quark index for quark */
  int combo_qk_index[MAX_QK][MAX_COMBO]; /* Quark index for combinations */
  int num_combo[MAX_QK];  /* Number of quarks to combine */
  Real combo_coeff[MAX_QK][MAX_COMBO]; /* Coefficients of linear combination */
  int naik_index[MAX_QK];            /* Naik term index for quark */
  quark_source snk_qs[MAX_QK];       /* Sink description for quark */
  int saveflag_q[MAX_QK];	     /* what to do for saving KS prop */
  char savefile_q[MAX_QK][MAXFILENAME];
  /* Mesons */
  int num_pair;                      /* Number of mesons */
  int qkpair[MAX_PAIR][2];           /* Indices of quarks in a meson */
  int do_meson_spect[MAX_PAIR];
  int saveflag_m[MAX_PAIR];          /* Save flag for meson correlator */
  char savefile_m[MAX_PAIR][MAXFILENAME]; /* File for meson correlator */
  int r_offset_m[MAX_PAIR][4];               /* Shift of origin for meson correlator */
  int num_corr_m[MAX_PAIR];                  /* Number of correlators for a meson */
  int num_corr_report[MAX_PAIR];           /* Number of correlators to report for a meson */
  char meson_label[MAX_PAIR][MAX_CORR][MAX_MESON_LABEL];
  char mom_label[MAX_PAIR][MAX_CORR][MAX_MOM_LABEL];
  int corr_index[MAX_PAIR][MAX_CORR];
  int spin_taste_snk[MAX_PAIR][MAX_CORR];
  int meson_phase[MAX_PAIR][MAX_CORR];
  Real meson_factor[MAX_PAIR][MAX_CORR];
  int corr_mom[MAX_PAIR][MAX_CORR][3];
  char corr_parity[MAX_PAIR][MAX_CORR][3];
  /* Baryons */
  int num_triplet;                   /* Number of baryons */
  int do_baryon_spect[MAX_TRIPLET];
  int num_corr_b[MAX_TRIPLET];             /* Number of correlators for a baryon */
  int saveflag_b[MAX_TRIPLET];            /* Save flag for baryon correlator */
  int baryon_type_snk[MAX_TRIPLET][MAX_CORR];
  char savefile_b[MAX_TRIPLET][MAXFILENAME]; /* File for baryon correlator */
  int qktriplet[MAX_TRIPLET][3];     /* Indices of quarks in a baryon */
  int r_offset_b[MAX_TRIPLET][4];          /* Shift of origin for baryon correlator */
  char baryon_label[MAX_TRIPLET][MAX_CORR][MAX_MESON_LABEL];
  int baryon_phase[MAX_TRIPLET][MAX_CORR];
  Real baryon_factor[MAX_TRIPLET][MAX_CORR];
  /* GB Baryons */
  #ifdef GB_BARYON
    int  num_oct;                               /* Number of ks quark octets */
    int  qk_oct[MAX_OCTET][8];                  /* List of ks quark octets */
    int  num_gb_triplet;                        /* Number of Golterman-Bailey baryons */
    int  qk8triplet[MAX_TRIPLET][3];            /* Nndices of octets used in gb baryon */
    int  qk8num_d[MAX_TRIPLET];                 /* Number of baryon d quarks */
    int  qk8num_s[MAX_TRIPLET];                 /* Number of baryon s quarks */
    int  do_gbbaryon_spect[MAX_TRIPLET];        /* Whether to do this correlator or not */
    int  gb_spintaste[MAX_TRIPLET];             /* Spin-taste phase construction at src/sink */
    int  r_offset_gb[MAX_TRIPLET][4];           /* Origin offset for correlator */
    int  snkmom_gb[MAX_TRIPLET][3];             /* Sink momentum */
    char snkpar_gb[MAX_TRIPLET][3];             /* Sink parity */
    int  num_corr_gb[MAX_TRIPLET];              /* Number of correlators for a gb baryon */
    char gbbaryon_label[MAX_TRIPLET][MAX_CORR][MAX_MESON_LABEL];
    int  gbbaryon_phase[MAX_TRIPLET][MAX_CORR];
    Real gbbaryon_factor[MAX_TRIPLET][MAX_CORR];
    int gbbaryon_src[MAX_TRIPLET][MAX_CORR]; /* gb op at source */
    int gbbaryon_snk[MAX_TRIPLET][MAX_CORR]; /* gb op at sink */
    int  gb_parity_flip_snk[MAX_TRIPLET][MAX_CORR]; /* Option to flip sink tie up epsilon location */
    int  saveflag_gb[MAX_TRIPLET];              /* Save flag for gb baryon correlator */
    char savefile_gb[MAX_TRIPLET][MAXFILENAME]; /* File for gb baryon correlator */
    char gb_mom_label[MAX_TRIPLET][MAX_CORR][MAX_MOM_LABEL]; /* gb momentum label */
    int  gb_corr_mom[MAX_TRIPLET][MAX_CORR][3];              /* gb correlator momentum */
    char gb_corr_parity[MAX_TRIPLET][MAX_CORR][3];           /* gb correlator parity */
    short gb_wall[MAX_TRIPLET][MAX_CORR]; /* Baryon with wall tieup or point tieup */
    short gb_corner[MAX_TRIPLET][MAX_CORR]; /* Baryon over corner or full cube */
  #endif
  /* Filenames */
  char startfile[MAXFILENAME];  /* Gauge file */
  char start_u1file[MAXFILENAME]; /* U(1) gauge file */
  char start_ape_file[MAXFILENAME]; /* APE links file */
  char savefile[MAXFILENAME];
  char save_u1file[MAXFILENAME]; /* U(1) gauge file */
  char save_ape_file[MAXFILENAME]; /* APE links file */
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
  char scratchstem_w[MAXFILENAME];
}  params;

//  short do_corr[MAX_PAIR][MAX_MESON][MAX_MESON_MOMENTUM];

#endif /* _PARAMS_H */
