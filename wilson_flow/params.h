#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */

/* structure for passing simulation parameters to each node */
typedef struct {
  int stopflag;   /* 1 if it is time to stop */
  
  /* INITIALIZATION PARAMETERS */
  int nx,ny,nz,nt;  /* lattice dimensions */

  /* FLOW PARAMETERS */
  Real stepsize; /* wilson flow time integration step size */
  Real stoptime; /* maximum flow time, -1 means auto-determined */
  int exp_order; /* where to end series expansion of exponential */
  char flow_description[20]; /* type of flow (wilson, symanzik) */
  int stapleflag; /* what type of action to use */

  int startflag;  /* what to do for beginning lattice */
  int saveflag;   /* what to do with lattice at end */
  char startfile[MAXFILENAME], savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
}  params;

#endif /* _PARAMS_H */
