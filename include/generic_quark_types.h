#ifndef _GENERIC_QUARK_TYPES_H
#define _GENERIC_QUARK_TYPES_H

#include "../include/macros.h"
#include "../include/precision.h"

/* Structure defining quark inversion parameters for most inverters */
typedef struct {
  int min;            /* minimum number of iterations */
  int max;            /* maximum number of iterations per restart */
  int nrestart;       /* maximum restarts */
  int parity;         /* EVEN, ODD, or EVENANDODD (for some inverters) */
  int start_flag;     /* 0: use a zero initial guess; 1: use dest */
  int nsrc;           /* Number of source vectors */
  Real resid;        /* desired residual - 
			 normalized as sqrt(r*r)/sqrt(src_e*src_e */
  Real size_r;       /* resulting residual */
  int converged;      /* returned 0 if not converged; 1 if converged */
  field_offset wv1;   /* ugly wilson_vector temporary */
  field_offset wv2;   /* ugly wilson_vector temporary */
  field_offset wv3;   /* ugly wilson_vector temporary */
  field_offset wv4;   /* ugly wilson_vector temporary */
                      /* Add further parameters as needed...  */
} quark_invert_control;

/* Structures required for specific inverters */

/* Structure defining parameters of Dirac matrix for clover inversion */
/* To be passed through to inverter. */
typedef struct {
  Real Kappa;        /* hopping */
  Real Clov_c;       /* Perturbative clover coeff */
  Real U0;           /* Tadpole correction to Clov_c */
  field_offset work_f_mn;       /* ugly su3_matrix temporary */
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
  int color;          /* source color required */
  int spin;           /* source spin required */
  int type;           /* source type for most source builders */
  char descrp[30];    /* alpha description for most */
  int wall_cutoff;    /* half size of box for w_source_h */
  int parity;         /* even or odd sites for w_source_h */
  Real r0;           /* source size for w_source */
  int x0,y0,z0,t0;    /* source coordinates for most */ 
  int src_pointer ;   /* smearing function (for the moment, only
		         clover_finite_p_vary/create_wilson_source.c) */
} wilson_quark_source;

#endif /* _GENERIC_QUARK_TYPES_H */


