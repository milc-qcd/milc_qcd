/********************  map_milc_to_qdp_Fd.c ***************************/
/* MILC Version 7 */

/* C-code include file for map_milc_to_qdp_F.c */
/* Defines precision-specific mappings between MILC and QDP for the
   case that the QDP precision is single and MILC precision is double */

/* 12/01/06 C. DeTar adapted from stuff_qdp.c */

#include <string.h>
#include <qdp.h>
#include "../include/generic_qdp.h"
#include "lattice.h"
#include "../include/macros.h"
#include "../include/complex.h"
#include "../include/su3.h"

#define make_set_qdp_from_site(T, TYPE, MILCTYPE) \
void \
set_F_##T##_from_site(QDP_F_##TYPE *dest, field_offset src, int parity) \
{ \
  int i; \
  site *s; \
  QLA_D_##TYPE t; \
  QLA_F_##TYPE *temp; \
  temp = QDP_F_expose_##T (dest); \
  FORSOMEPARITY(i,s,parity) { \
    memcpy((void *)&t, F_PT(s,src), sizeof(QLA_D_##TYPE)); \
    QLA_FD_##T##_eq_##T (&temp[i], &t); \
  } \
  QDP_F_reset_##T (dest); \
}

#define make_set_site_from_qdp(T, TYPE, MILCTYPE) \
void \
set_site_from_F_##T(field_offset dest, QDP_F_##TYPE *src, int parity) \
{ \
  int i; \
  site *s; \
  QLA_D_##TYPE t; \
  QLA_F_##TYPE *temp; \
  temp = QDP_F_expose_##T (src); \
  FORSOMEPARITY(i,s,parity) { \
    QLA_DF_##T##_eq_##T (&t, &temp[i]); \
    memcpy((void *)F_PT(s,dest), (void *)&t, sizeof(QLA_D_##TYPE)); \
  } \
  QDP_F_reset_##T (src); \
}

#define make_set_qdp_from_field(T, TYPE, MILCTYPE) \
void \
set_F_##T##_from_field(QDP_F_##TYPE *dest, MILCTYPE *src, int parity)	\
{ \
  int i; \
  site *s; \
  QLA_D_##TYPE t; \
  QLA_F_##TYPE *temp; \
  temp = QDP_F_expose_##T (dest); \
  FORSOMEPARITY(i,s,parity) { \
    memcpy((void *)&t, (void *)&src[i], sizeof(QLA_D_##TYPE)); \
    QLA_FD_##T##_eq_##T (&temp[i],&t); \
  } \
  QDP_F_reset_##T (dest); \
}

#define make_set4_qdp_from_site(T, TYPE, MILCTYPE) \
void \
set4_F_##T##_from_site(QDP_F_##TYPE *dest[], field_offset src, int parity)\
{ \
  int i, dir; \
  site *s; \
  QLA_D_##TYPE t; \
  QLA_F_##TYPE *temp; \
  for(dir=0; dir<4; dir++) { \
    temp = QDP_F_expose_##T (dest[dir]); \
    FORSOMEPARITY(i,s,parity) { \
      memcpy((void *)&t, (char *)F_PT(s,src)+dir*sizeof(MILCTYPE), sizeof(QLA_D_##TYPE)); \
    QLA_FD_##T##_eq_##T (&temp[i],&t); \
    } \
    QDP_F_reset_##T (dest[dir]); \
  } \
}

#define make_set_field_from_qdp(T, TYPE, MILCTYPE) \
void \
 set_field_from_F_##T(MILCTYPE *dest, QDP_F_##TYPE *src, int parity) \
{ \
  int i; \
  site *s; \
  QLA_D_##TYPE t; \
  QLA_F_##TYPE *temp; \
  temp = QDP_F_expose_##T (src); \
  FORSOMEPARITY(i,s,parity) { \
    QLA_DF_##T##_eq_##T (&t, &temp[i]); \
    memcpy((QLA_D_##TYPE *)&dest[i], (void *)&t, sizeof(QLA_D_##TYPE)); \
  } \
  QDP_F_reset_##T (src); \
}

#define make_set4_qdp_from_field(T, TYPE, MILCTYPE) \
void \
 set4_F_##T##_from_field(QDP_F_##TYPE *dest[], MILCTYPE *src, int parity) \
{ \
  int i, dir; \
  site *s; \
  QLA_D_##TYPE t; \
  QLA_F_##TYPE *temp; \
  for(dir=0; dir<4; dir++) { \
    temp = QDP_F_expose_##T (dest[dir]); \
    FORSOMEPARITY(i,s,parity) { \
      memcpy((void *)&t, (void *)&src[4*i+dir], sizeof(QLA_D_##TYPE)); \
      QLA_FD_##T##_eq_##T (&temp[i],&t); \
    } \
    QDP_F_reset_##T (dest[dir]); \
  } \
}

#define make_set4_field_from_qdp(T, TYPE, MILCTYPE) \
void \
set4_field_from_F_##T(MILCTYPE *dest, QDP_F_##TYPE *src[], int parity)	\
{ \
  int i, dir; \
  site *s; \
  QLA_D_##TYPE t; \
  QLA_F_##TYPE *temp; \
  for(dir=0; dir<4; dir++) { \
    temp = QDP_F_expose_##T (src[dir]); \
    FORSOMEPARITY(i,s,parity) { \
      QLA_DF_##T##_eq_##T (&t, &temp[i]); \
      memcpy((void *)&dest[4*i+dir], (void *)&t, sizeof(QLA_D_##TYPE)); \
    } \
    QDP_F_reset_##T (src[dir]); \
  } \
}

copy_types(make_set_qdp_from_site);

copy_types(make_set_site_from_qdp);

copy_types(make_set_qdp_from_field);

copy_types(make_set_field_from_qdp);

copy_types(make_set4_field_from_qdp);

copy_types(make_set4_qdp_from_field);

copy_types(make_set4_qdp_from_site);

/* map_milc_to_qdp_Fd.c */
