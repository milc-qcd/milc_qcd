#ifndef _PARAMS_H
#define _PARAMS_H

#include "defines.h"
#include "../include/macros.h"  /* For MAXFILENAME */
#include "../include/generic_quark_types.h"
#include "../include/generic_ks.h" /* For quark_source */
#include "../include/generic_wilson.h"  /* For quark_source */
#include "../include/gammatypes.h"

#define MAX_MASS_PBP 8
#define MAX_SOURCE 8
#define MAX_SET 8
#define MAX_PROP 32
#define MAX_QK 64
#define MAX_PAIR 512
#define MAX_TRIPLET 32
#define MAX_QKPAIR_LABEL 32
#define MAX_MESON 32
#define MAX_SPECTRUM_REQUEST 512
#define MAX_MESON_LABEL 32
#define MAX_BARYON_LABEL 32
#define MAX_MESON_MOMENTUM 100
#define MAX_MOM_LABEL 16
#define MAX_CORR 200
#define STATIC_TYPE 0
#define KS_TYPE 1
#define PROP_TYPE 0
#define QUARK_TYPE 1
#define BASE_SOURCE_PARENT -1

enum checktype { CHECK_NO,  CHECK_YES, CHECK_SOURCE_ONLY };

/* structure for passing simulation parameters to each node */
typedef struct {
  int stopflag;   /* 1 if it is time to stop */
  /* INITIALIZATION PARAMETERS */
  int nx,ny,nz,nt;	/* lattice dimensions */
#ifdef FIX_NODE_GEOM
  int node_geometry[4];  /* Specifies fixed "nsquares" (i.e. 4D
			    hypercubes) for the compute nodes in each
			    coordinate direction.  Must be divisors of
			    the lattice dimension */
#ifdef FIX_IONODE_GEOM
  int ionode_geometry[4]; /* Specifies fixed "nsquares" for I/O
			     partitions in each coordinate direction,
			     one I/O node for each square.  The I/O
			     node is at the origin of the square.
			     Must be divisors of the node_geometry. */
#endif
#endif
  int iseed;
  char job_id[MAXFILENAME]; /* Usually encoded by scripts */

  /*  REPEATING BLOCK */
  int startflag;	/* what to do for beginning lattice */
  Real u0;
  int coord_origin[4];  /* Origin of coordinates for KS phases and time_bc */
  int fixflag;    /* whether to gauge fix */
  int saveflag;	/* what to do for saving lattice */
  Real staple_weight;
  int ape_iter;
  int num_pbp_masses;   /* Number of masses for pbp calculation */
  quark_invert_control qic_pbp[MAX_MASS_PBP];
  int prec_pbp;         /* Precision of the pbp calculation (1 or 2) */
  int npbp_reps;     /* Number of random sources for pbp calculation */
  ks_param ksp_pbp[MAX_MASS_PBP];
  int num_base_source;  /* Number of base sources */
  quark_source base_src_qs[MAX_SOURCE];
  int num_modified_source;       /* Number of modified sources */
  quark_source_sink_op src_qs_op[MAX_SOURCE];
  int parent_source[MAX_SOURCE];      /* base_source or source index */
  int num_set;  /* number of sets */
  int source[MAX_SET];      /* index of modified source for this set */
  int num_prop[MAX_SET]; /* number of propagators in a set */
  quark_source src_qs[MAX_SET];
  int prop_type[MAX_PROP]; /* 0 static 1 KS */
  int set[MAX_PROP];  /* The set to which the propagator belongs */
  int begin_prop[MAX_SET]; /* index of beginning propagator in this set */
  int end_prop[MAX_SET]; /* index of ending propagator in this set */
  int startflag_ks[MAX_PROP];	/* what to do for beginning KS prop */
  int saveflag_ks[MAX_PROP];	/* what to do for saving KS prop */
  char mass_label[MAX_PROP][32];
  quark_source_sink_op snk_qs_op[MAX_QK];
  ks_param ksp[MAX_PROP];
  quark_invert_control qic[MAX_PROP];
  int check[MAX_PROP];         /* True -> run the inverter */
  Real bdry_phase[MAX_PROP][4];      /* For twisted boundary conditions */
  int num_qk;	                     /* number of quarks */
  int parent_type[MAX_QK];           /* propagator type: quark or propagator */
  int prop_for_qk[MAX_QK];           /* Propagator or quark index for quark */
  int naik_index[MAX_QK];            /* Naik term index for quark */
  quark_source snk_qs[MAX_QK];       /* Sink description for quark */
  int saveflag_q[MAX_QK];	     /* what to do for saving KS prop */
  int num_pair;                      /* Number of mesons */
  int qkpair[MAX_PAIR][2];           /* Indices of quarks in a meson */
  int num_triplet;                   /* Number of baryons */
  int qktriplet[MAX_TRIPLET][3];     /* Indices of quarks in a baryon */
  int do_meson_spect[MAX_PAIR];      
  int do_baryon_spect[MAX_TRIPLET];
  int saveflag_m[MAX_PAIR];          /* Save flag for meson correlator */
  char savefile_m[MAX_PAIR][MAXFILENAME]; /* File for meson correlator */
  int saveflag_b[MAX_TRIPLET];            /* Save flag for baryon correlator */
  char savefile_b[MAX_TRIPLET][MAXFILENAME]; /* File for baryon correlator */
  int r_offset_m[MAX_PAIR][4];               /* Shift of origin for meson correlator */
  int r_offset_b[MAX_TRIPLET][4];          /* Shift of origin for baryon correlator */
  int num_corr_m[MAX_PAIR];                  /* Number of correlators for a meson */
  int num_corr_b[MAX_TRIPLET];             /* Number of correlators for a baryon */
  int num_corr_report[MAX_PAIR];           /* Number of correlators to report for a meson */
  char meson_label[MAX_PAIR][MAX_CORR][MAX_MESON_LABEL];
  char mom_label[MAX_PAIR][MAX_CORR][MAX_MOM_LABEL];
  int corr_index[MAX_PAIR][MAX_CORR];
  int spin_taste_snk[MAX_PAIR][MAX_CORR];
  int baryon_type_snk[MAX_TRIPLET][MAX_CORR];
  int meson_phase[MAX_PAIR][MAX_CORR];
  Real meson_factor[MAX_PAIR][MAX_CORR];
  char baryon_label[MAX_TRIPLET][MAX_CORR][MAX_MESON_LABEL];
  int baryon_phase[MAX_TRIPLET][MAX_CORR];
  Real baryon_factor[MAX_TRIPLET][MAX_CORR];
  int corr_mom[MAX_PAIR][MAX_CORR][3];
  char corr_parity[MAX_PAIR][MAX_CORR][3];
  char startfile[MAXFILENAME];  /* Gauge file */
  char savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
  char scratchstem_w[MAXFILENAME];
  char startfile_ks[MAX_PROP][MAXFILENAME];
  char savefile_ks[MAX_PROP][MAXFILENAME];
  char savefile_q[MAX_QK][MAXFILENAME];
}  params;

//  short do_corr[MAX_PAIR][MAX_MESON][MAX_MESON_MOMENTUM];

#endif /* _PARAMS_H */
