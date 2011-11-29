#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "../include/generic_quark_types.h"
#include "../include/generic_ks.h" /* For quark_source */
#include "../include/generic_wilson.h"  /* For quark_source */
#include "../include/gammatypes.h"

#define MAX_SOURCE 24
#define MAX_PROP 16
#define MAX_QK 64
#define MAX_PAIR 128
#define MAX_QKPAIR_LABEL 32
#define MAX_MESON 32
#define MAX_SPECTRUM_REQUEST 512
#define MAX_MESON_LABEL 32
#define MAX_MESON_MOMENTUM 100
#define MAX_MOM_LABEL 16
#define MAX_CORR 200
#define CLOVER_TYPE 0
#define KS_TYPE 1
#define KS4_TYPE 2
#define PROP_TYPE 0
#define QUARK_TYPE 1
#define BASE_SOURCE_PARENT -1

enum checktype { CHECK_NO, CHECK_YES, CHECK_SOURCE_ONLY };

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
  int fixflag;    /* whether to gauge fix */
  int saveflag;	/* what to do for saving lattice */
  Real staple_weight;
  int ape_iter;
  Real u0;
  int coord_origin[4];  /* Origin of coordinates for KS phases and time_bc */
  quark_invert_control qic[MAX_PROP];
  int num_base_source;       /* Number of modified sources */
  quark_source base_src_qs[MAX_SOURCE];
  int num_modified_source;       /* Number of sources modifications */
  quark_source_sink_op src_qs_op[MAX_SOURCE];
  int parent_source[MAX_SOURCE];      /* base_source or source index */
  int num_prop; /* number of propagators */
  int prop_type[MAX_PROP]; /* 0 clover 1 KS */
  int source[MAX_PROP]; /* source index for this propagator */
  int startflag_w[MAX_PROP];	/* what to do for beginning wilson prop */
  int saveflag_w[MAX_PROP];	/* what to do for saving wilson prop */
  quark_source src_qs[MAX_PROP];
  dirac_clover_param dcp[MAX_PROP];
  char kappa_label[MAX_PROP][32];
  int startflag_ks[MAX_PROP];	/* what to do for beginning wilson prop */
  int saveflag_ks[MAX_PROP];	/* what to do for saving KS prop */
  ks_param ksp[MAX_PROP];
  char mass_label[MAX_PROP][32];
  int check[MAX_PROP];
  Real bdry_phase[MAX_PROP][4];      /* For twisted boundary conditions */
  int num_qk;	              /* number of quarks */
  int parent_type[MAX_QK];      /* propagator type: quark or propagator */
  int prop_for_qk[MAX_QK];           /* Propagator or quark index for quark */
  int naik_index[MAX_QK];            /* Naik term index for quark */
  int saveflag_q[MAX_QK];	/* what to do for saving wilson prop */
  quark_source_sink_op snk_qs_op[MAX_QK];
  int num_pair;
  int qkpair[MAX_PAIR][2];
  int do_meson_spect[MAX_PAIR];
  int do_baryon_spect[MAX_PAIR];
  int do_closed_hadron[MAX_PAIR];
  int do_open_meson[MAX_PAIR];
  int saveflag_c[MAX_PAIR];
  char savefile_c[MAX_PAIR][MAXFILENAME];
  int t_offset[MAX_PAIR];
  int r_offset[MAX_PAIR][3];
  int num_corr[MAX_PAIR];
  int num_corr_report[MAX_PAIR];
  char meson_label[MAX_PAIR][MAX_CORR][MAX_MESON_LABEL];
  char mom_label[MAX_PAIR][MAX_CORR][MAX_MOM_LABEL];
  int corr_index[MAX_PAIR][MAX_CORR];
  int gam_src[MAX_PAIR][MAX_CORR], gam_snk[MAX_PAIR][MAX_CORR];
  int corr_phase[MAX_PAIR][MAX_CORR];
  Real corr_factor[MAX_PAIR][MAX_CORR];
  int corr_mom[MAX_PAIR][MAX_CORR][3];
  char corr_parity[MAX_PAIR][MAX_CORR][3];
  char startfile[MAXFILENAME];
  char savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
  char scratchstem_w[MAXFILENAME];
  char startfile_w[MAX_PROP][MAXFILENAME];
  char savefile_w[MAX_PROP][MAXFILENAME];
  char startfile_ks[MAX_PROP][MAXFILENAME];
  char savefile_ks[MAX_PROP][MAXFILENAME];
  char savefile_q[MAX_QK][MAXFILENAME];
}  params;

//  short do_corr[MAX_PAIR][MAX_MESON][MAX_MESON_MOMENTUM];

#endif /* _PARAMS_H */
