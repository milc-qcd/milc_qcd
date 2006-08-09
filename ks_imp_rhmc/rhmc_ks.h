#ifndef _RHMC_KS_H
#define _RHMC_KS_H
/************************ rhmc_ks.h **********************************
*									*
*  Macros and declarations for RHMC algorithm routines                  *
*/


int ks_ratinv(	/* Return value is number of iterations taken */
    field_offset src,   /* source vector (type su3_vector) */
    su3_vector **psim,  /* solution vectors */
    Real mass,          /* quark mass */
    Real *roots,        /* the roots */
    int order,          /* order of rational function approx */
    int niter,          /* maximal number of CG interations */
    Real rsqmin,        /* desired residue squared */
    int parity,         /* parity to be worked on */
    Real *final_rsq_ptr /* final residue squared */
    );

int ks_rateval(
    su3_vector *dest,   /* answer vector */
    field_offset src,   /* source vector (for a_0 term) */
    su3_vector **psim,  /* solution vectors  from multiroot CG */
    Real *residues,     /* the residues */
    int order,          /* order of approximation */
    int parity          /* parity to be worked on */
    );

#endif /* _RHMC_KS_H */
