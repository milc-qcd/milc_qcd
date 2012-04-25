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

/* Inverter types */
enum cl_cg_t { BICG, CG, MR, HOP };

int congrad_cl(int niter,Real rsqmin,Real *final_rsq_ptr);

int bicgilu_cl_site(     /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int bicgilu_cl_field_cpu(    /* Return value is number of iterations taken */
    wilson_vector *src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int bicgilu_cl_field_gpu(    /* Return value is number of iterations taken */
    wilson_vector *src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

#ifdef USE_CL_GPU
#define bicgilu_cl_field bicgilu_cl_field_gpu
#else
#define bicgilu_cl_field bicgilu_cl_field_cpu
#endif

int bicgilu_w_site(       /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int bicgilu_w_field(     /* Return value is number of iterations taken */
    wilson_vector *src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int cgilu_cl_site(       /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int cgilu_cl_field(       /* Return value is number of iterations taken */
    wilson_vector *src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int cgilu_w_site(        /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int cgilu_w_field(       /* Return value is number of iterations taken */
    wilson_vector *src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int hopilu_cl_site(      /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int hopilu_cl_field(     /* Return value is number of iterations taken */
    wilson_vector *src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int hopilu_w_site(       /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int hopilu_w_field(      /* Return value is number of iterations taken */
    wilson_vector *src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int mrilu_cl_site(       /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int mrilu_cl_field(      /* Return value is number of iterations taken */
    wilson_vector *src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int mrilu_w_site(        /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int mrilu_w_field(       /* Return value is number of iterations taken */
    wilson_vector *src,  /* type wilson_vector (source vector - OVERWRITTEN!)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

/* cl_solver_utilities.c */
Real relative_residue(wilson_vector *p, wilson_vector *q, int parity);

#include "../include/comdefs.h"

double ilu_xfm_source(
     wilson_vector *t_dest,
     wilson_vector *r,
     wilson_vector *my_mp,
     Real Kappa,
     int *is_startede,
     msg_tag *tage[]
     );

void ilu_DRD(
  wilson_vector *t_dest,
  wilson_vector *my_mp,
  wilson_vector *tmp,
  wilson_vector *tmpo,
  int isign,
  msg_tag* tago[],
  int *is_startedo,
  msg_tag* tage[],
  int *is_startede
  );

void ilu_xfm_dest(
     wilson_vector *t_dest,
     wilson_vector *my_mp,
     Real Kappa,
     int *is_startedo,
     msg_tag *tago[]);

void map_dwp_to_dcp(dirac_clover_param *dcp, dirac_wilson_param *dwp);

/* f_mu_nu.c */

void f_mu_nu(su3_matrix f_mn[],int mu,int nu);

/* make_clov2.c */

typedef struct { complex tr[2][15]; } triangular;
typedef struct { Real di[2][6]; } diagonal;
typedef struct {
  triangular *clov;
  diagonal *clov_diag;
  triangular *clov_raw;
  diagonal *clov_diag_raw;
  double trlogA;
  Real Clov_c;
  int valid_clov;
} clover;

void mult_this_ldu_site(
  clover *my_clov,
  field_offset src,   /* type wilson_vector RECAST AS wilson_block_vector */
  field_offset dest,  /* type wilson_vector RECAST AS wilson_block_vector */
  int parity
			);


/* make_clov.c routines for any clover term */
clover *create_clov(void);
void compute_clov(clover *my_clov, Real Clov_c);
void compute_clovinv(clover *my_clov, int parity);
void mult_this_ldu_field(
  clover *my_clov,
  wilson_vector *src,  /* RECAST as wilson_block_vector */
  wilson_vector *dest, /* RECAST as wilson_block_vector */
  int parity
  );
void tr_sigma_this_ldu_mu_nu_site( clover *my_clov, field_offset mat, int mu, 
				   int nu );
void invalidate_this_clov(clover *my_clov);
void free_this_clov(clover *my_clov);
void destroy_this_clov(clover **my_clov);

/* make_clov.c routines for single clover term */
void make_clov(Real Clov_c);
double make_clovinv(int parity);
void free_clov(void);
void mult_ldu_site(
  field_offset src,   /* type wilson_vector RECAST AS wilson_block_vector */
  field_offset dest,  /* type wilson_vector RECAST AS wilson_block_vector */
  int parity );
void mult_ldu_field(
  wilson_vector *src,
  wilson_vector *dest,
  int parity
  );
void invalidate_clov(void);
#endif /* _GENERIC_CLOVER_H */
