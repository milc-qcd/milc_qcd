#ifndef _GENERIC_U1
#define _GENERIC_U1

/************************ generic_u1.h **********************************
*									*
*  Macros and declarations for generic_u1 routines                      *
*  This header is for codes that call generic_u1 routines               *
*  MIMD version 7 							*
*									*/

#include "../include/complex.h"

/*****************************************************************************/
/* Everything needed to build the U(1)-modified fermion links for a set of
   U(1) fields.  The caller owns the fields themselves: in use case (a)
   they come from generate_u1_gauge_field(), in use case (b) they are
   read in.  Either way they are plain vector potentials, indexed
   u1_A[iU1][4*i+dir], as create_u1_A_field() allocates them. */

typedef struct {
  int nU1;           /* Number of U(1) vector potentials */
  Real **u1_A;       /* [nU1] fields, each 4*sites_on_node Reals */
  int n_naiks;       /* Naik epsilon table, as for create_fermion_links_from_site */
  double *eps_naik;
  int i_naik;        /* Which Naik epsilon fn_mass corresponds to */
} u1_field_info_t;

/* u1avlink.c */
void u1avlink(double *sLink, double *tLink, double charge);

/* u1link.c */
Real *create_u1_A_field(void);
Real **create_u1_A_field_array(int nU1_A);
void destroy_u1_A_field(Real *A);
void destroy_u1_A_field_array(int nU1_A, Real **A);
void u1phase_on(Real charge, Real *A);
void u1phase_off(void);

/* u1plaq.c */
complex u1ploop(Real charge);
void u1plaq(Real *ssplq,Real *stplq,Real charge);

/* u1_gauge_field.c */

/* Choice of gauge for the generated vector potential */
typedef enum {
  U1_FEYNMAN = 0,   /* No gauge fixing */
  U1_LORENZ  = 1,   /* Project out the longitudinal mode k_mu       */
  U1_COULOMB = 2    /* Project out the longitudinal spatial mode    */
} u1_gauge_fix_t;

/* Treatment of the modes with vanishing spatial momentum */
typedef enum {
  U1_QED_NONE = 0,  /* No special treatment                                */
  U1_QED_L    = 1,  /* Drop all components when vec(k) = 0                 */
  U1_QED_TL   = 2   /* Drop the time component when vec(k) = 0 and, if     */
                    /*  k_t = 0 as well, the spatial components also       */
} u1_zero_mode_t;

/* Choice of generation algorithm.  This one stays a compile-time macro.
   An external caller that does not define GAUGE_METHOD gets the general
   algorithm, which is the only one that honors the gauge and scheme
   arguments below. */
#define U1_GENERAL          0  /* Gauge and scheme selected at run time */
#define U1_MILC_ORIGINAL    1  /* Legacy MILC code.  Known to be incorrect */
#define U1_HATTON_REVISION  2  /* Dan Hatton's revision */

#ifndef GAUGE_METHOD
#define GAUGE_METHOD U1_GENERAL
#endif

void initialize_u1_gauge_field(u1_gauge_fix_t gauge, u1_zero_mode_t scheme);
void generate_u1_gauge_field(Real *u1_A);
void generate_u1_gauge_field_array(int nU1_A, Real *u1_A_array[]);
void finalize_u1_gauge_field(void);
void show_u1_gauge_field_opts(u1_gauge_fix_t gauge, u1_zero_mode_t scheme);

#endif /* _GENERIC_U1 */
