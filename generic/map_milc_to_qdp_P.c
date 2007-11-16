/********************  map_milc_to_qdp_P.c ***************************/
/* MILC Version 7 */

/* C-code include file for map_milc_to_qdp_F.c and map_milc_to_qdp_D.c */
/* Defines precision-specific mappings between MILC and QDP for the
   case that the precisions of both are equal  */

/* 12/01/06 C. DeTar adapted from stuff_qdp.c */

#include <qdp.h>
#include "../include/generic_qdp.h"
#include <string.h>
#include "lattice.h"
#include "../include/macros.h"
#include "../include/complex.h"
#include "../include/su3.h"

#define make_set_qdp_from_site(P, T, TYPE, MILCTYPE) \
void \
set_##P##_##T##_from_site(QDP_##TYPE *dest, field_offset src, int parity) \
{ \
  int i; \
  site *s; \
  QLA_##TYPE *temp; \
  temp = QDP_expose_##T (dest); \
  FORSOMEPARITY(i,s,parity) {   \
    memcpy((void *)&temp[i], F_PT(s,src), sizeof(QLA_##TYPE)); \
  } \
  QDP_reset_##T (dest); \
}

#define make_set_site_from_qdp(P, T, TYPE, MILCTYPE) \
void \
set_site_from_##P##_##T(field_offset dest, QDP_##TYPE *src, int parity)	\
{ \
  int i; \
  site *s; \
  QLA_##TYPE *temp; \
  temp = QDP_expose_##T (src); \
  FORSOMEPARITY(i,s,parity) {	\
    memcpy(F_PT(s,dest), (void *)&temp[i], sizeof(QLA_##TYPE)); \
  } \
  QDP_reset_##T (src); \
}

#define make_set_qdp_from_field(P, T, TYPE, MILCTYPE) \
void \
set_##P##_##T##_from_field(QDP_##TYPE *dest, MILCTYPE *src, int parity) \
{ \
  int i; \
  site *s; \
  QLA_##TYPE *temp; \
  temp = QDP_expose_##T (dest); \
  FORSOMEPARITY(i,s,parity) {  \
    memcpy((void *)&temp[i], (void *)&src[i], sizeof(QLA_##TYPE)); \
  } \
  QDP_reset_##T (dest); \
}

#define make_set4_qdp_from_site(P, T, TYPE, MILCTYPE)	\
void \
set4_##P##_##T##_from_site(QDP_##TYPE *dest[], field_offset src, int parity) \
{ \
  int i, dir; \
  site *s; \
  QLA_##TYPE *temp; \
  for(dir=0; dir<4; dir++) { \
    temp = QDP_expose_##T (dest[dir]); \
    FORSOMEPARITY(i,s,parity) {						\
      memcpy((void *)&temp[i], (char *)F_PT(s,src)+dir*sizeof(MILCTYPE), sizeof(QLA_##TYPE)); \
    } \
    QDP_reset_##T (dest[dir]); \
  } \
}

#define make_set_field_from_qdp(P, T, TYPE, MILCTYPE) \
void \
set_field_from_##P##_##T(MILCTYPE *dest, QDP_##TYPE *src, int parity)	\
{ \
  int i; \
  site *s; \
  QLA_##TYPE *temp; \
  temp = QDP_expose_##T (src); \
  FORSOMEPARITY(i,s,parity) {   \
    memcpy((void *)&dest[i], (void *)&temp[i], sizeof(QLA_##TYPE)); \
  } \
  QDP_reset_##T (src); \
}

#define make_set4_qdp_from_field(P, T, TYPE, MILCTYPE) \
void \
set4_##P##_##T##_from_field(QDP_##TYPE *dest[], MILCTYPE *src, int parity) \
{ \
  int i, dir; \
  site *s; \
  QLA_##TYPE *temp; \
  for(dir=0; dir<4; dir++) { \
    temp = QDP_expose_##T (dest[dir]); \
    FORSOMEPARITY(i,s,parity) {	\
      memcpy((void *)&temp[i], (void *)&src[4*i+dir], sizeof(QLA_##TYPE)); \
    } \
    QDP_reset_##T (dest[dir]); \
  } \
}

#define make_set4_field_from_qdp(P, T, TYPE, MILCTYPE) \
void \
set4_field_from_##P##_##T(MILCTYPE *dest, QDP_##TYPE *src[], int parity) \
{ \
  int i, dir; \
  site *s; \
  QLA_##TYPE *temp; \
  for(dir=0; dir<4; dir++) { \
    temp = QDP_expose_##T (src[dir]); \
    FORSOMEPARITY(i,s,parity) { \
      memcpy((void *)&dest[4*i+dir], (void *)&temp[i], sizeof(QLA_##TYPE)); \
    } \
    QDP_reset_##T (src[dir]); \
  } \
}

copy_types(make_set_qdp_from_site);

copy_types(make_set_site_from_qdp);

copy_types(make_set_qdp_from_field);

copy_types(make_set_field_from_qdp);

copy_types(make_set4_field_from_qdp);

copy_types(make_set4_qdp_from_field);

copy_types(make_set4_qdp_from_site);

/* map_milc_to_qdp_P.c */
