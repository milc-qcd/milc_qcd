/********************  map_milc_to_qdp_F.c ***************************/
/* MILC Version 7 */

/* Defines routines for precision-specific mappings between MILC and
   QDP for the case that the QDP precision is double and the
   prevailing MILC precision is single or double */

/* We assume that the MILC precision for site and field variables
   is uniformly set by PRECISION, but the QDP precision can be set
   in each call if so desired */

/* 12/01/06 C. DeTar adapted from stuff_to_qdp.c */

#undef QDP_Precision
#define QDP_Precision 'D'

#if (PRECISION == 2 )

#define copy_types(macro) \
macro(D, R, Real, Real); \
macro(D, V, ColorVector, su3_vector); \
macro(D, H, HalfFermion, half_wilson_vector); \
macro(D, D, DiracFermion, wilson_vector); \
macro(D, M, ColorMatrix, su3_matrix);

#include "map_milc_to_qdp_P.c"

#else

#define copy_types(macro) \
macro(R, Real, Real); \
macro(V, ColorVector, su3_vector); \
macro(H, HalfFermion, half_wilson_vector); \
macro(D, DiracFermion, wilson_vector); \
macro(M, ColorMatrix, su3_matrix);

#include "map_milc_to_qdp_Df.c"

#endif

/* map_milc_to_qdp_F.c */
