#ifndef QOP_MILC_H
#define QOP_MILC_H

#include "../include/su3.h"

/* For MILC test implementation of SciDAC Level 3 routines */

#if QOP_Precision == 1

#define QOP_ColorVector_struct QOP_F3_ColorVector_struct
#define QOP_DiracFermion_struct QOP_F3_DiracFermion_struct
#define QOP_GaugeField_struct QOP_F3_GaugeField_struct
#define QOP_Force_struct QOP_F3_Force_struct
#define QOP_FermionLinksAsqtad_struct QOP_F3_FermionLinksAsqtad_struct

#else

#define QOP_ColorVector_struct QOP_D3_ColorVector_struct
#define QOP_DiracFermion_struct QOP_D3_DiracFermion_struct
#define QOP_GaugeField_struct QOP_D3_GaugeField_struct
#define QOP_Force_struct QOP_D3_Force_struct
#define QOP_FermionLinksAsqtad_struct QOP_D3_FermionLinksAsqtad_struct

#endif

struct QOP_ColorVector_struct {
  su3_vector *v;
  int evenodd;
};

struct QOP_DiracFermion_struct {
  wilson_vector *d;
  int evenodd;
};

struct QOP_GaugeField_struct {
  su3_matrix *g;
  int evenodd;
};

struct QOP_Force_struct {
  su3_matrix *f;
  int evenodd;
};

struct QOP_FermionLinksAsqtad_struct {
  struct QOP_GaugeField_struct *fat;
  struct QOP_GaugeField_struct *lng;
  int evenodd;
};

#include <qop.h>

#endif /* QOP_MILC_H */
