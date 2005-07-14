#ifndef _GENERIC_CLOVER_H
#define _GENERIC_CLOVER_H
/************************ generic_clover.h ******************************
*									*
*  Macros and declarations for generic_clover routines                  *
*  This header is for codes that call generic_clover routines           *
*  MIMD version 7 							*
*									*
*/

#include "../include/generic_quark_types.h"
#include "../include/macros.h"

int cgilu_cl(            /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int hopilu_cl(           /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int bicgilu_cl(          /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int mrilu_cl(            /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int congrad_cl(int niter,Real rsqmin,Real *final_rsq_ptr);

void f_mu_nu(field_offset f_mn,int mu,int nu);

void make_clov(Real Clov_c,field_offset f_mn);
double make_clovinv(int parity);
void tr_sigma_ldu_mu_nu( field_offset mat, int mu, int nu );
void free_clov();
void mult_ldu(
  field_offset src,   /* type wilson_vector RECAST AS wilson_block_vector */
  field_offset dest,  /* type wilson_vector RECAST AS wilson_block_vector */
  int parity );
void mult_ldu_on_temp(
  wilson_vector *src,
  wilson_vector *dest,
  int parity
  );
#endif /* _GENERIC_CLOVER_H */
