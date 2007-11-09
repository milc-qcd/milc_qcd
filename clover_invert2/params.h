#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "../include/generic_quark_types.h" /* For wilson_quark_source */
#include "../include/generic_wilson.h"
#include "../include/gammatypes.h"

#define MAX_QK 6
#define MAX_PAIR 8
#define MAX_MESON 32
#define MAX_SPECTRUM_REQUEST 512
#define MAX_MESON_LABEL 32
#define MAX_SINK_LABEL 32
#define MAX_MESON_MOMENTUM 20
#define MAX_MOM_LABEL 16

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
  quark_invert_control qic;
  int num_qk;	/* number of quarks */
  int startflag_w[MAX_QK];	/* what to do for beginning wilson vector */
  int saveflag_w[MAX_QK];	/* what to do for saving wilson vector */
  dirac_clover_param dcp[MAX_QK];
  wilson_quark_source src_wqs[MAX_QK];
  int check[MAX_QK];
  Real d1[MAX_QK];
  int num_pair;
  int qkpair[MAX_PAIR][2];
  wilson_quark_source snk_wqs[MAX_QK];
  char snk_label[MAX_QK][MAX_SINK_LABEL];
  int do_point_meson_spect[MAX_PAIR];
  int do_smear_meson_spect[MAX_PAIR];
  int do_rot_meson_spect[MAX_PAIR];
  int do_baryon_spect[MAX_PAIR];
  int saveflag_c[MAX_PAIR];
  char savefile_c[MAX_PAIR][MAXFILENAME];
  int num_meson_report;
  char meson_label[MAX_MESON][MAX_MESON_LABEL];
  int num_meson;
  int meson_index[MAX_MESON];
  int gam_src[MAX_MESON], gam_snk[MAX_MESON];
  complex meson_phase[MAX_MESON];
  int num_mom;
  int meson_mom[MAX_MESON_MOMENTUM][3];
  int num_mom_report;
  char mom_label[MAX_MESON_MOMENTUM][MAX_MOM_LABEL];
  int mom_index[MAX_MESON_MOMENTUM];
  char startfile[MAXFILENAME];
  char savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
  char scratchstem_w[MAXFILENAME];
  char startfile_w[MAX_QK][MAXFILENAME];
  char savefile_w[MAX_QK][MAXFILENAME];
}  params;


#endif /* _PARAMS_H */
