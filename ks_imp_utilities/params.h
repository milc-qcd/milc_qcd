#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "defines.h"

/* structure for passing simulation parameters to each node */
typedef struct {
  int stopflag;   /* 1 if it is time to stop */
  /* INITIALIZATION PARAMETERS */
  int nx,ny,nz,nt;  /* lattice dimensions */
  int iseed;	/* for random numbers */
  /*  REPEATING BLOCK */
  Real mass; /* gauge coupling, quark masses */
  Real u0; /* tadpole parameter */
  int niter; 	/* maximum number of c.g. iterations */
  int nrestart;	/* maximum number of c.g. restarts */
  Real rsqprop;  /* for deciding on convergence */
  int startflag;  /* what to do for beginning lattice */
  int saveflag;   /* what to do with lattice at end */
  int savelongflag, savefatflag;  /* same for longlinks and fatlinks */
  int srcflag; /* what to do for source lattice */
  int ansflag; /* what to do for answer lattice */
  
  char startfile[MAXFILENAME],savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
  char savelongfile[MAXFILENAME],savefatfile[MAXFILENAME];
  char srcfile[MAXFILENAME],ansfile[MAXFILENAME];
  int inverttype;
}  params;

#endif /* _PARAMS_H */
