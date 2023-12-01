#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */

/* structure for passing simulation parameters to each node */
typedef struct {
  int stopflag;   /* 1 if it is time to stop */

  /* INITIALIZATION PARAMETERS */
  int nx,ny,nz,nt;  /* lattice dimensions */

  /* FLOW PARAMETERS */
#ifdef ANISOTROPY
  Real ani;
#endif
  Real stepsize; /* wilson flow time integration step size */
  Real stoptime; /* maximum flow time, -1 means auto-determined */
  int exp_order; /* where to end series expansion of exponential */
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  Real local_tol; /* local tolerance for adaptive integrators */
#endif

#ifdef REGIONS
  Real stepsize_bulk; /* wilson flow time integration step size for bulk*/
  Real stoptime_bulk; /* maximum flow time for bulk, -1 means auto-determined */
  int exp_order_bulk; /* where to end series expansion of exponential */
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  Real local_tol_bulk; /* local tolerance for adaptive integrators */
#endif

  Real stepsize_bdry; /* wilson flow time integration step size for boundary*/
  Real stoptime_bdry; /* maximum flow time for boundary, -1 means auto-determined */
  int exp_order_bdry; /* where to end series expansion of exponential */
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  Real local_tol_bdry; /* local tolerance for adaptive integrators */
#endif
#endif

#ifdef SPHALERON
  Real qthr_bulk;
  int minqthr_bulk;
  int maxnflow_bulk;

  Real qthr_bdry;
  int minqthr_bdry;
  int maxnflow_bdry;

  Real qs_tol;
#endif
#ifdef BLOCKING
  Real block_1to2_time;
  Real block_2to4_time;
#endif

  char flow_description[20]; /* type of flow (wilson, symanzik, zeuthen) */
  int stapleflag; /* what type of action to use */

#ifdef REGIONS
  char flow_description_bulk[20]; /* type of flow (wilson, symanzik, zeuthen) */
  int stapleflag_bulk; /* what type of action to use */

  char flow_description_bdry[20]; /* type of flow (wilson, symanzik, zeuthen) */
  int stapleflag_bdry; /* what type of action to use */
#endif

  int startflag;  /* what to do for beginning lattice */
  int saveflag;   /* what to do with lattice at end */
  char startfile[MAXFILENAME], savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
}  params;

#endif /* _PARAMS_H */
