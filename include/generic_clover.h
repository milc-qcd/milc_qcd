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

void f_mu_nu(su3_matrix f_mn[],int mu,int nu);

typedef struct { complex tr[2][15]; } triangular;
typedef struct { Real di[2][6]; } diagonal;
typedef struct {
  triangular *clov;
  diagonal *clov_diag;
} clover;

/* make_clov.c routines for any clover term */
clover *create_clov(void);
void compute_clov(clover *my_clov, Real Clov_c);
double compute_clovinv(clover *my_clov, int parity);
void mult_this_ldu_field(
  clover *my_clov,
  wilson_vector *src,  /* RECAST as wilson_block_vector */
  wilson_vector *dest, /* RECAST as wilson_block_vector */
  int parity
  );
void free_this_clov(clover *my_clov);

/* make_clov.c routines for single clover term */
void make_clov(Real Clov_c);
double make_clovinv(int parity);
void tr_sigma_ldu_mu_nu_site( field_offset mat, int mu, int nu );
void free_clov();
void mult_ldu_site(
  field_offset src,   /* type wilson_vector RECAST AS wilson_block_vector */
  field_offset dest,  /* type wilson_vector RECAST AS wilson_block_vector */
  int parity );
void mult_ldu_field(
  wilson_vector *src,
  wilson_vector *dest,
  int parity
  );
#endif /* _GENERIC_CLOVER_H */
