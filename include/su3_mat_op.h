/* Temporary include file.
   These routined should eventually be included in su3.h
   Started: 071005
   Last modified: 081121 */


#ifndef SU3_MAT_OP_H_
#define SU3_MAT_OP_H_

/* MILC includes */
#include "../include/complex.h"
#include "../include/su3.h"

#if (PRECISION==1)
#define U3_UNIT_ANALYTIC_EPS 1.0e-6
#define U3_ROOT_INV_NORM_EPS 1.0e-6
#else
#define U3_UNIT_ANALYTIC_EPS 1.0e-14
#define U3_ROOT_INV_NORM_EPS 1.0e-9
#endif

#define U3_UNIT_DER_EPS 1.0e-6

// If one runs into small eigenvalues when calculating
// the HISQ fermion force, defining this option allows to
// regularize them and preven large spikes in the force

// #define HISQ_FORCE_FILTER 5.0e-5 // Now defined in project Make_template

#define U3_ROOT_INV_MAX_ITER 100
#define MILC_AB_PI 3.14159265358979323846264338328
#define MILC_AB_TPI 6.28318530717958647692528676656

#define U3_UNIT_RAT_NTERMS 14
#ifndef SU3_MAT_OP_NO_STORAGE
static Real c_l_U3_UNIT_RAT[ U3_UNIT_RAT_NTERMS+1 ] = {
         0.0850910,
         2.413330975,
         6.257184884e-01,
         2.925707925e-01,
         1.737405612e-01,
         1.166359792e-01,
         8.372555094e-02,
         6.216038074e-02,
         4.652496186e-02,
         3.423610040e-02,
         2.404754621e-02,
         1.545550091e-02,
         8.436481876e-03,
         3.419245947e-03,
         1.138166539e-03  };
static Real d_l_U3_UNIT_RAT[ U3_UNIT_RAT_NTERMS+1 ] = {
         0.0,
         1.361747338e+01,
         3.135687028e+00,
         1.213113539e+00,
         5.596349298e-01,
         2.752627333e-01,
         1.364115846e-01,
         6.543005714e-02,
         2.923946484e-02,
         1.164228894e-02,
         3.887745892e-03,
         9.937321442e-04,
         1.684882417e-04,
         1.585925699e-05,
         5.914114023e-07  };
#endif


#endif /* SU3_MAT_OP_H_ */
