#ifndef _GF_GLOBALS_H
#define _GF_GLOBALS_H
/****************************** gf_globals.h ********************************/
/* global parameters for gradient flow */

/* Maximum number of stages for 2N-storage RK schemes */
#define MAX_RK_2N 5

/* Integrator parameters */
/* 2N-storage schemes */
#if INTEGRATOR==INTEGRATOR_LUSCHER || INTEGRATOR==INTEGRATOR_CK
// number of stages
EXTERN int N_2N;
// A, B coefficients
EXTERN double A_2N[MAX_RK_2N];
EXTERN double B_2N[MAX_RK_2N];
#endif


#endif /* _GF_GLOBALS_H */
