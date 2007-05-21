#ifndef _GENERIC_QUARK_TYPES_H
#define _GENERIC_QUARK_TYPES_H

#include "../include/macros.h"
#include "../include/precision.h"

/* Structure defining quark inversion parameters for most inverters */
typedef struct {
  int prec;           /* precision of the inversion 1 = single; 2 = double */
  int min;            /* minimum number of iterations (being phased out) */
  int max;            /* maximum number of iterations per restart */
  int nrestart;       /* maximum restarts */
  int parity;         /* EVEN, ODD, or EVENANDODD (for some inverters) */
  int start_flag;     /* 0: use a zero initial guess; 1: use dest */
  int nsrc;           /* Number of source vectors */
  Real resid;         /* desired residual - 
			 normalized as sqrt(r*r)/sqrt(src_e*src_e) */
  Real relresid;      /* desired relative residual */
  Real final_rsq;     /* Final true (absolute) residual */
  Real final_relrsq;  /* Final relative residual */
  Real size_r;        /* resulting cumulative residual */
  Real size_relr;     /* resulting cumulative relative residual */
  int converged;      /* returned 0 if not converged; 1 if converged */
                      /* Add further parameters as needed...  */
} quark_invert_control;

/* Structures required for specific inverters */

/* Structure defining parameters of Dirac matrix for clover inversion */
/* To be passed through to inverter. */
typedef struct {
  Real Kappa;        /* hopping */
  Real Clov_c;       /* Perturbative clover coeff */
  Real U0;           /* Tadpole correction to Clov_c */
} dirac_clover_param;

/* Same for Wilson case */
typedef struct {
  Real Kappa;        /* hopping */
} dirac_wilson_param;

/* Same for plain KS case */
typedef struct {
  Real mass;
} dirac_ks_param;

/* Structure defining Wilson (or clover) quark source */
/* There must be a color and spin member */
/* Add other members to suit the generic_wilson code 
   that builds the source.  Ignore the members you don't need. */
typedef struct {
  int type;           /* source type for most source builders */
  char descrp[30];    /* alpha description for most */
  int color;          /* source color */
  int spin;           /* source spin  */
  int wall_cutoff;    /* half size of box for w_source_h */
  int parity;         /* even or odd sites for w_source_h */
  Real r0;            /* source size for gaussian, width for gauge invt  */
  int iters;          /* iterations for gauge invariant source */
  int x0,y0,z0,t0;    /* source coordinates for most */ 
  char source_file[MAXFILENAME]; /* file name for some sources */
  int src_pointer ;   /* smearing function (for the moment, only
		         clover_finite_p_vary/create_wilson_source.c) */
} wilson_quark_source;

#endif /* _GENERIC_QUARK_TYPES_H */


