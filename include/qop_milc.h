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

struct QOP_ColorVector_struct {
  fsu3_vector *v;
  int evenodd;
};

struct QOP_DiracFermion_struct {
  fwilson_vector *d;
  int evenodd;
};

struct QOP_GaugeField_struct {
  fsu3_matrix *g;
  int evenodd;
};

struct QOP_Force_struct {
  fsu3_matrix *f;
  int evenodd;
};

#else

#define QOP_ColorVector_struct QOP_D3_ColorVector_struct
#define QOP_DiracFermion_struct QOP_D3_DiracFermion_struct
#define QOP_GaugeField_struct QOP_D3_GaugeField_struct
#define QOP_Force_struct QOP_D3_Force_struct
#define QOP_FermionLinksAsqtad_struct QOP_D3_FermionLinksAsqtad_struct

struct QOP_ColorVector_struct {
  dsu3_vector *v;
  int evenodd;
};

struct QOP_DiracFermion_struct {
  dwilson_vector *d;
  int evenodd;
};

struct QOP_GaugeField_struct {
  dsu3_matrix *g;
  int evenodd;
};

struct QOP_Force_struct {
  dsu3_matrix *f;
  int evenodd;
};

#endif

struct QOP_FermionLinksAsqtad_struct {
  struct QOP_GaugeField_struct *fat;
  struct QOP_GaugeField_struct *lng;
  int evenodd;
};

/* Precision conversion routines for unpacking the data members above 
   into flat arrays */
/* qop_milc_utilities.c */

su3_matrix * create_links_from_qop_milc_F(fsu3_matrix *src);
void destroy_links_from_qop_milc_F(su3_matrix *g);
su3_matrix *create_links_from_qop_milc_D(dsu3_matrix *src);
void destroy_links_from_qop_milc_D(su3_matrix *g);
su3_vector *create_latvec_from_qop_milc_F(fsu3_vector *src);
void destroy_latvec_from_qop_milc_F(su3_vector *v);
su3_vector *create_latvec_from_qop_milc_D(dsu3_vector *src);
void destroy_latvec_from_qop_milc_D(su3_vector *v);
void copy_latvec_to_qop_milc_F( fsu3_vector *dest, su3_vector *src);
void copy_latvec_to_qop_milc_D( dsu3_vector *dest, su3_vector *src);

#include <qop.h>

#endif /* QOP_MILC_H */
