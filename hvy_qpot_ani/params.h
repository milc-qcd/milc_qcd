#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */

/* structure for passing simulation parameters to each node */
typedef struct {
  int stopflag;	/* 1 if it is time to stop */
  /* INITIALIZATION PARAMETERS */
  int nx,ny,nz,nt;  /* lattice dimensions */
  Real u0;
  int fixflag; /* whether to gauge fix */

  /*  details of the correlations to be considered; 
   *  the time direction for the correlation function does not 
   *  necessarily agree with the time direction TUP of the lattice */
  short cor_dir;
  int min_ct;
  int max_x, max_y, max_z, max_t;
#ifdef ANISOTROPY
  short ani_dir; /* direction of anisotropy */
  Real ani_xiq; /* bare quark anisotropy */
#ifdef ONEDIM_ANISO_TEST
  Real iso_xiq; /* bare quark isotropic link factor for debugging */
#endif
#endif
  Real max_r;
  int off_axis_flag;	/* off-axis Wilson loops or not? */

  /*  details of the smearing to be used */
  int no_smear_level;	/* number of smearing levels (<=5) */
  int smear_num[5];	/* the number of smearing iterations */
#ifdef HYP_SMEARING /* BOTH in 3D and 4D HYP smearing */
  Real hyp_alpha1;      /* parameters for 4D HYP smearing */
  Real hyp_alpha2;
  Real hyp_alpha3;
#else /* APE smearing */
  int ape_iter;
  Real staple_weight;
#endif
  Real smear_fac;       /* smearing factor = weight of direct link */
#if (defined APE_1D_SMEARING || defined APE_1D2_SMEARING)
  int stap_dir;
#endif
#ifdef NEW_HVY_POT
int hqp_alg;
#endif

  /*  details of starting and saving */
  int startflag;	/* what to do for beginning lattice */
  int saveflag;	/* what to do with lattice at end */
  char startfile[MAXFILENAME],savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
}  params;

#endif /* _PARAMS_H */
