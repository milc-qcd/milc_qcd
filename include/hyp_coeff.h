#ifndef _HYP_COEFF_H
#define _HYP_COEFF_H

#include "../include/complex.h"

typedef struct {
  Real alpha1;
  Real alpha2;
  Real alpha3;
  int proj_method;
  int hits;
} hyp_coeffs_t;

// project to U(3) in HYP action
#define HYP_U3 0
// project to SU(3) in HYP action by projecting to U(3)
//   and then removing the phase
#define HYP_SU3 1
// project to SU(3) in HYP action by trace maximization
#define HYP_SU3_TR_MAX 2


// set HYP coefficients
void set_hyp_coeff( hyp_coeffs_t *hc, Real a1, Real a2, Real a3 );

// set group projection method: U(3), SU(3), SU(3) by Tr Max,
// last parameter is the number of hits for the third method,
// will be ignored by first two methods
void set_hyp_proj_method( hyp_coeffs_t *hc, int method, int hits );

// HYP smearing
// if dir_exclude==NODIR -- usual 4D HYP smearing
// if dir_exclude==XUP or YUP or ZUP or TUP, 3D smearing
void load_hyp_links(su3_matrix *U_link, su3_matrix *hyp_link, int dir_exclude, hyp_coeffs_t *hc);


#endif /* _HYP_COEFF_H */

