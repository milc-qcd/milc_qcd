#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "../include/generic_quark_types.h"
#include "../include/generic_ks.h" /* For ks_quark_source */
#include "../include/generic_wilson.h"  /* For wilson_quark_source */
#include "../include/gammatypes.h"

#define MAX_QK 6
#define MAX_PAIR 8
#define MAX_QKPAIR_LABEL 32
#define MAX_MESON 32
#define MAX_SPECTRUM_REQUEST 512
#define MAX_MESON_LABEL 32
#define MAX_SINK_LABEL 32
#define MAX_MESON_MOMENTUM 100
#define MAX_MOM_LABEL 16
#define MAX_CORR 200
#define MAX_CORR_LABEL 
#define CLOVER_TYPE 0
#define KS_TYPE 1

/* structure for passing simulation parameters to each node */
typedef struct {
  int stopflag;   /* 1 if it is time to stop */
  /* INITIALIZATION PARAMETERS */
  int nx,ny,nz,nt;	/* lattice dimensions */
  char job_id[MAXFILENAME]; /* Usually encoded by scripts */

  /*  REPEATING BLOCK */
  int startflag;	/* what to do for beginning lattice */
  int fixflag;    /* whether to gauge fix */
  int saveflag;	/* what to do for saving lattice */
  int scratchflag;
  Real u0;
  quark_invert_control qic;
  int num_qk;	/* number of quarks */
  int qk_type[MAX_QK];          /* 0 clover 1 KS */
  int startflag_w[MAX_QK];	/* what to do for beginning wilson vector */
  int saveflag_w[MAX_QK];	/* what to do for saving wilson vector */
  dirac_clover_param dcp[MAX_QK];
  wilson_quark_source src_wqs[MAX_QK];
  int startflag_ks[MAX_QK];	/* what to do for beginning wilson vector */
  int saveflag_ks[MAX_QK];	/* what to do for saving wilson vector */
  ks_param ksp[MAX_QK];
  ks_quark_source src_ksqs[MAX_QK];
  int check[MAX_QK];
  Real d1[MAX_QK];
  int num_pair;
  int qkpair[MAX_PAIR][2];
  wilson_quark_source snk_wqs[MAX_QK];
  int do_point_meson_spect[MAX_PAIR];
  int do_smear_meson_spect[MAX_PAIR];
  int do_rot_meson_spect[MAX_PAIR];
  int do_baryon_spect[MAX_PAIR];
  int saveflag_c[MAX_PAIR];
  char savefile_c[MAX_PAIR][MAXFILENAME];
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
  char startfile_w[MAX_QK][MAXFILENAME];
  char savefile_w[MAX_QK][MAXFILENAME];
  char startfile_ks[MAX_QK][MAXFILENAME];
  char savefile_ks[MAX_QK][MAXFILENAME];
}  params;

//  short do_corr[MAX_PAIR][MAX_MESON][MAX_MESON_MOMENTUM];

#endif /* _PARAMS_H */
