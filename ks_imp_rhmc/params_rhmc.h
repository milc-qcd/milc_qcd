#ifndef _PARAMS_RHMC_H
#define _PARAMS_RHMC_H

#define NMASS 4

typedef struct {
  int  y[NMASS];
  int  z[NMASS];
  Real m[NMASS];
  int order;
  Real *res;
  Real *pole;
} params_ratfunc;

typedef struct {
  params_ratfunc MD;
  params_ratfunc GR;
  params_ratfunc FA;
  Real naik_term_epsilon;
} params_rhmc;
#endif /* _PARAMS_RHMC_H */
