#include <string.h>
#include <qdp.h>
#include "lattice.h"
#include "../include/macros.h"
#include "../include/complex.h"
#include "../include/su3.h"

#define make_set_qdp_from_site(T, TYPE, MILCTYPE) \
void \
set_##T##_from_site(QDP_##TYPE *dest, field_offset src) \
{ \
  int i; \
  site *s; \
  QLA_##TYPE *temp; \
  temp = QDP_expose_##T (dest); \
  FORALLSITES(i,s) { \
    memcpy((void *)&temp[i], F_PT(s,src), sizeof(QLA_##TYPE)); \
  } \
  QDP_reset_##T (dest); \
}

#define make_set_site_from_qdp(T, TYPE, MILCTYPE) \
void \
set_site_from_##T(field_offset dest, QDP_##TYPE *src) \
{ \
  int i; \
  site *s; \
  QLA_##TYPE *temp; \
  temp = QDP_expose_##T (src); \
  FORALLSITES(i,s) { \
    memcpy(F_PT(s,dest), (void *)&temp[i], sizeof(QLA_##TYPE)); \
  } \
  QDP_reset_##T (src); \
}

#define make_set_qdp_from_field(T, TYPE, MILCTYPE) \
void \
set_##T##_from_field(QDP_##TYPE *dest, MILCTYPE *src) \
{ \
  int i; \
  site *s; \
  QLA_##TYPE *temp; \
  temp = QDP_expose_##T (dest); \
  FORALLSITES(i,s) { \
    memcpy((void *)&temp[i], (void *)&src[i], sizeof(QLA_##TYPE)); \
  } \
  QDP_reset_##T (dest); \
}

#define make_set4_qdp_from_site(T, TYPE, MILCTYPE) \
void \
set4_##T##_from_site(QDP_##TYPE *dest[], field_offset src) \
{ \
  int i, dir; \
  site *s; \
  QLA_##TYPE *temp; \
  for(dir=0; dir<4; dir++) { \
    temp = QDP_expose_##T (dest[dir]); \
    FORALLSITES(i,s) { \
      memcpy((void *)&temp[i], (void *)F_PT(s,src)+dir*sizeof(MILCTYPE), sizeof(QLA_##TYPE)); \
    } \
    QDP_reset_##T (dest[dir]); \
  } \
}

#define make_set_field_from_qdp(T, TYPE, MILCTYPE) \
void \
set_field_from_##T(MILCTYPE *dest, QDP_##TYPE *src) \
{ \
  int i; \
  site *s; \
  QLA_##TYPE *temp; \
  temp = QDP_expose_##T (src); \
  FORALLSITES(i,s) { \
    memcpy((void *)&dest[i], (void *)&temp[i], sizeof(QLA_##TYPE)); \
  } \
  QDP_reset_##T (src); \
}

#define make_set4_qdp_from_field(T, TYPE, MILCTYPE) \
void \
set4_##T##_from_field(QDP_##TYPE *dest[], MILCTYPE *src) \
{ \
  int i, dir; \
  site *s; \
  QLA_##TYPE *temp; \
  for(dir=0; dir<4; dir++) { \
    temp = QDP_expose_##T (dest[dir]); \
    FORALLSITES(i,s) { \
      memcpy((void *)&temp[i], (void *)&src[4*i+dir], sizeof(QLA_##TYPE)); \
    } \
    QDP_reset_##T (dest[dir]); \
  } \
}

#define make_set4_field_from_qdp(T, TYPE, MILCTYPE) \
void \
set4_field_from_##T(MILCTYPE *dest, QDP_##TYPE *src[]) \
{ \
  int i, dir; \
  site *s; \
  QLA_##TYPE *temp; \
  for(dir=0; dir<4; dir++) { \
    temp = QDP_expose_##T (src[dir]); \
    FORALLSITES(i,s) { \
      memcpy((void *)&dest[4*i+dir], (void *)&temp[i], sizeof(QLA_##TYPE)); \
    } \
    QDP_reset_##T (src[dir]); \
  } \
}

#define copy_types(macro) \
macro(I, Int, int); \
macro(R, Real, Real); \
macro(V, ColorVector, su3_vector); \
macro(H, HalfFermion, half_wilson_vector); \
macro(D, DiracFermion, wilson_vector); \
macro(M, ColorMatrix, su3_matrix);

copy_types(make_set_qdp_from_site);

copy_types(make_set_site_from_qdp);

copy_types(make_set_qdp_from_field);

copy_types(make_set_field_from_qdp);

copy_types(make_set4_field_from_qdp);

copy_types(make_set4_qdp_from_field);

copy_types(make_set4_qdp_from_site);

