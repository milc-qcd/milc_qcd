#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */

/* structure for passing simulation parameters to each node */
typedef struct {
  int stopflag;	/* 1 if it is time to stop */
  /* INITIALIZATION PARAMETERS */
  int nx,ny,nz,nt;  /* lattice dimensions */
  /*  REPEATING BLOCK */
  int max_t;
  int max_x;
  int ape_iter;
  Real u0;
  Real staple_weight;
  int no_smear_level;	/* number of smearing levels (<=5) */
  int smear_num[5];	/* the number of smearing iterations */
  int off_axis_flag;	/* off-axis Wilson loops or not? */
  int startflag;	/* what to do for beginning lattice */
  int saveflag;	/* what to do with lattice at end */
  Real beta;	/* gauge coupling */
  Real smear_fac;	/* smearing factor = weight of direct link */
  char startfile[MAXFILENAME],savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
}  params;

#endif /* _PARAMS_H */
