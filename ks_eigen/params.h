#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/imp_ferm_links.h"

/* structure for passing simulation parameters to each node */
typedef struct {
  int stopflag;   /* 1 if it is time to stop */
  /* INITIALIZATION PARAMETERS */
  int nx,ny,nz,nt;  /* lattice dimensions */
  int iseed;	/* for random numbers */
  /*  REPEATING BLOCK */
  Real mass; /* gauge coupling, quark mass */
  Real u0; /* tadpole parameter */
  Real staple_weight;
  int ape_iter;
  int niter; 	/* maximum number of c.g. iterations */
  int nrestart; 	/* maximum number of c.g. restarts */
  Real rsqmin,rsqprop;  /* for deciding on convergence */
  ks_eigen_param eigen_param; /* Parameters for eigensolver */
  int ks_eigen_saveflag; /* eigenvector file type */
  char ks_eigen_savefile[MAXFILENAME]; /* eigenvector output file name */
  
  int startflag;  /* what to do for beginning lattice */
  char startfile[MAXFILENAME]; /* starting lattice file name */
}  params;

#endif /* _PARAMS_H */
