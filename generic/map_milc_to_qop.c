OBSOLETE
/*********************** map_milc_to_qop.c **************************/
/* Functions for mapping MILC data layouts to raw QOP layouts       */
/* C. DeTar 10/19/2005                                              */

#include "generic_includes.h"
#include "../include/generic_qop.h"

/* Create empty raw links */

#define make_create_raw4(P, T, MILCTYPE) \
MILCTYPE ** \
create_raw4_##P##_##T (void){ \
  MILCTYPE **raw = NULL; \
  int dir; \
  raw = (MILCTYPE **)malloc(4*sizeof(MILCTYPE *)); \
  FORALLUPDIR(dir){ \
    raw[dir] = (MILCTYPE *)malloc(sites_on_node*sizeof(MILCTYPE)); \
    if(raw[dir] == NULL){ \
      printf("create4_raw: No room\n"); \
      return NULL; \
    } \
  } \
  return raw; \
}

/* Destroy raw links */

#define make_destroy_raw4(P, T, MILCTYPE) \
void \
destroy_raw4_##P##_##T (MILCTYPE *raw[]){ \
  int dir; \
  FORALLUPDIR(dir){ \
    if(raw[dir] != NULL) \
      free(raw[dir]); \
  } \
  free(raw); \
}

/* Create empty raw field */

#define make_create_raw(P, T, MILCTYPE) \
MILCTYPE * \
create_raw_##P##_##T(void){ \
  MILCTYPE *raw = NULL; \
  raw = (MILCTYPE *)malloc(sites_on_node*sizeof(MILCTYPE)); \
  if(raw == NULL){ \
    printf("create_raw: No room\n"); \
    return NULL; \
  } \
  return raw; \
}

/* Destroy raw field */

#define make_destroy_raw(P, T, MILCTYPE) \
void \
destroy_raw_##P##_##T (MILCTYPE *raw){ \
  if(raw != NULL) free(raw); \
}

/* Copy types with possible conversion */

/* Convert (or copy) MILC types between specific and prevailing precision */

#if (PRECISION==1)

static void 
f2p_mat(su3_matrix *dest, fsu3_matrix *src){
  memcpy((void *)dest, (void *)src, sizeof(fsu3_matrix));
}

static void 
p2f_mat(fsu3_matrix *dest, su3_matrix *src){
  memcpy((void *)dest, (void *)src, sizeof(fsu3_matrix));
}

static void 
d2p_mat(su3_matrix *dest, dsu3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
p2d_mat(dsu3_matrix *dest, su3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
f2p_vec(su3_vector *dest, fsu3_vector *src){
  memcpy((void *)dest, (void *)src, sizeof(fsu3_vector));
}

static void 
p2f_vec(fsu3_vector *dest, su3_vector *src){
  memcpy((void *)dest, (void *)src, sizeof(fsu3_vector));
}

static void 
d2p_vec(su3_vector *dest, dsu3_vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

static void 
p2d_vec(dsu3_vector *dest, su3_vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

#else

static void 
f2p_mat(su3_matrix *dest, fsu3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

/* Convert (or copy) su3_matrix from prevailing to single precision */
static void 
p2f_mat(fsu3_matrix *dest, su3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
d2p_mat(su3_matrix *dest, dsu3_matrix *src){
  memcpy((void *)dest, (void *)src, sizeof(dsu3_matrix));
}

static void 
p2d_mat(dsu3_matrix *dest, su3_matrix *src){
  memcpy((void *)dest, (void *)src, sizeof(dsu3_matrix));
}

static void 
f2p_vec(su3_vector *dest, fsu3_vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

/* Convert (or copy) su3_vector from prevailing to single precision */
static void 
p2f_vec(fsu3_vector *dest, su3_vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

static void 
d2p_vec(su3_vector *dest, dsu3_vector *src){
  memcpy((void *)dest, (void *)src, sizeof(dsu3_vector));
}

static void 
p2d_vec(dsu3_vector *dest, su3_vector *src){
  memcpy((void *)dest, (void *)src, sizeof(dsu3_vector));
}

#endif

/* Conversions from prevailing MILC formats to specified formats */

#define copy_milc_to_F_G(d,s) p2f_mat(d,s);
#define copy_milc_to_D_G(d,s) p2d_mat(d,s);

static void
copy_milc_to_F_F(fsu3_matrix *dest, anti_hermitmat *src){
  su3_matrix t;
  uncompress_anti_hermitian( src, &t );
  p2f_mat( dest, &t );
}

static void
copy_milc_to_D_F(dsu3_matrix *dest, anti_hermitmat *src){
  su3_matrix t;
  uncompress_anti_hermitian( src, &t );
  p2d_mat( dest, &t );
}

#define copy_milc_to_F_V(d,s) p2f_vec(d,s);
#define copy_milc_to_D_V(d,s) p2d_vec(d,s);

/* Conversions from specified formats to prevailing MILC formats */

#define copy_F_G_to_milc(d,s) f2p_mat(d,s);
#define copy_D_G_to_milc(d,s) d2p_mat(d,s);

static void
copy_F_F_to_milc(anti_hermitmat *dest, fsu3_matrix *src){
  su3_matrix t;
  f2p_mat(&t, src);
  make_anti_hermitian( &t, dest );
}

static void
copy_D_F_to_milc(anti_hermitmat *dest, dsu3_matrix *src){
  su3_matrix t;
  d2p_mat(&t, src);
  make_anti_hermitian( &t, dest );
}

#define copy_F_V_to_milc(d,s) f2p_vec(d,s);
#define copy_D_V_to_milc(d,s) d2p_vec(d,s);

static void site_coords(int coords[4],site *s){
  coords[0] = s->x;
  coords[1] = s->y;
  coords[2] = s->z;
  coords[3] = s->t;
}

/* Map MILC site links to raw order */

#define make_create_raw4_from_site(P, T, MILC_RAWTYPE, MILC_SRCTYPE) \
MILC_RAWTYPE ** \
create_raw4_##P##_##T##_from_site(field_offset src, int milc_parity){ \
  int coords[4]; \
  int i,j,dir; \
  site *s; \
  MILC_RAWTYPE **raw; \
  MILC_SRCTYPE *tmp; \
  raw = create_raw4_##P##_##T (); \
  if(raw == NULL)return NULL; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("create_raw4_from_site: incompatible layout\n"); \
      return NULL; \
    } \
    j = QOP_node_index_raw_##T(coords, milc2qop_parity(milc_parity)); \
    FORALLUPDIR(dir){ \
      tmp = (MILC_SRCTYPE *)F_PT(s, src); \
      copy_milc_to_##P##_##T(raw[dir] + j, &tmp[dir]); \
    } \
  } \
  return raw; \
}

/* Map MILC links in site-major order to raw */

#define make_create_raw4_from_field(P, T, MILC_RAWTYPE, MILC_SRCTYPE) \
MILC_RAWTYPE ** \
create_raw4_##P##_##T##_from_field(MILC_SRCTYPE *src, int milc_parity){ \
  int coords[4]; \
  int i,j,dir; \
  site *s; \
  MILC_RAWTYPE **raw = NULL; \
  MILC_SRCTYPE *tmp; \
  raw = create_raw4_##P##_##T (); \
  if(raw == NULL)return NULL; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("create_raw4_from_field: incompatible layout\n"); \
      return NULL; \
    } \
    j = QOP_node_index_raw_##T (coords, milc2qop_parity(milc_parity)); \
    FORALLUPDIR(dir){ \
      tmp = src + 4*i; \
      copy_milc_to_##P##_##T(raw[dir] + j, &tmp[dir]); \
    } \
  } \
  return raw; \
}

/* Map MILC site field to raw */

#define make_create_raw_from_site(P, T, MILC_RAWTYPE, MILC_SRCTYPE) \
MILC_RAWTYPE * \
create_raw_##P##_##T##_from_site(field_offset src, int milc_parity){ \
  int coords[4]; \
  int i,j,dir; \
  site *s; \
  MILC_RAWTYPE *raw; \
  MILC_SRCTYPE *tmp; \
  raw = create_raw_##P##_##T(); \
  if(raw == NULL)return NULL; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("create_raw_from_site: incompatible layout\n"); \
      return NULL; \
    } \
    j = QOP_node_index_raw_##T (coords, milc2qop_parity(milc_parity)); \
    FORALLUPDIR(dir){ \
      tmp = (MILC_SRCTYPE *)F_PT(s, src); \
      copy_milc_to_##P##_##T(raw + j, tmp); \
    } \
  } \
  return raw; \
}

/* Map MILC field to raw */

#define make_create_raw_from_field(P, T, MILC_RAWTYPE, MILC_SRCTYPE) \
MILC_RAWTYPE * \
create_raw_##P##_##T##_from_field(MILC_SRCTYPE *src, int milc_parity){ \
  int coords[4]; \
  int i,j,dir; \
  site *s; \
  MILC_RAWTYPE *raw; \
  raw = create_raw_##P##_##T(); \
  if(raw == NULL)return NULL; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("create_raw_from_field: incompatible layout\n"); \
      return NULL; \
    } \
    j = QOP_node_index_raw_##T (coords, milc2qop_parity(milc_parity)); \
    FORALLUPDIR(dir){ \
      copy_milc_to_##P##_##T(raw + j, src + i); \
    } \
  } \
  return raw; \
}

/* Map raw links to MILC site structure */

#define make_unload_raw4_to_site(P, T, MILC_DSTTYPE, MILC_RAWTYPE) \
void \
unload_raw4_##P##_##T##_to_site(field_offset dest, MILC_RAWTYPE *raw[], \
         int milc_parity){ \
  int coords[4]; \
  int i,j,dir; \
  site *s; \
  MILC_DSTTYPE *tmp; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("unload_raw4_to_site: incompatible layout\n"); \
      terminate(1); \
    } \
    j = QOP_node_index_raw_##T(coords, milc2qop_parity(milc_parity)); \
    FORALLUPDIR(dir){ \
      tmp = (MILC_DSTTYPE *)F_PT(s, dest); \
      copy_##P##_##T##_to_milc(&tmp[dir], raw[dir] + j); \
    } \
  } \
}

#define make_unload_raw4_to_field(P, T, MILC_DSTTYPE, MILC_RAWTYPE) \
void \
unload_raw4_##P##_##T##_to_field(MILC_DSTTYPE *dest, MILC_RAWTYPE *raw[], \
         int milc_parity){ \
  int coords[4]; \
  int i,j,dir; \
  site *s; \
  MILC_DSTTYPE *tmp; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("unload_raw4_to_field: incompatible layout\n"); \
      terminate(1); \
    } \
    j = QOP_node_index_raw_##T(coords, milc2qop_parity(milc_parity)); \
    FORALLUPDIR(dir){ \
      tmp = dest + 4*i; \
      copy_##P##_##T##_to_milc(&tmp[dir], raw[dir] + j); \
    } \
  } \
}

#define make_unload_raw_to_site(P, T, MILC_DSTTYPE, MILC_RAWTYPE) \
void \
unload_raw_##P##_##T##_to_site(field_offset dest, MILC_RAWTYPE *raw, \
       int milc_parity){ \
  int coords[4]; \
  int i,j; \
  site *s; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("unload_raw_to_site: incompatible layout\n"); \
      terminate(1); \
    } \
    j = QOP_node_index_raw_##T(coords, milc2qop_parity(milc_parity)); \
    copy_##P##_##T##_to_milc((MILC_DSTTYPE *)F_PT(s,dest), raw + j); \
  } \
}

#define make_unload_raw_to_field(P, T, MILC_DSTTYPE, MILC_RAWTYPE) \
void \
unload_raw_##P##_##T##_to_field(MILC_DSTTYPE *dest, MILC_RAWTYPE *raw, \
       int milc_parity){ \
  int coords[4]; \
  int i,j; \
  site *s; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("unload_raw_V_to_field: incompatible layout\n"); \
      terminate(1); \
    } \
    j = QOP_node_index_raw_##T(coords, milc2qop_parity(milc_parity)); \
    copy_##P##_##T##_to_milc(dest + i, raw + j); \
  } \
}

#define make_load_links_and_mom_site(P, MILCTYPE, MILCFLOAT) \
void \
load_##P##_links_and_mom_site(QOP_##P##3_GaugeField **links, \
   QOP_##P##3_Force **mom, MILCTYPE ***rawlinks, MILCTYPE ***rawmom) \
{ \
  *rawlinks = create_raw4_##P##_G_from_site(F_OFFSET(link), EVENANDODD); \
  if(*rawlinks == NULL)terminate(1); \
  *links = QOP_##P##3_create_G_from_raw((MILCFLOAT **)(*rawlinks),QOP_EVENODD); \
  *rawmom = create_raw4_##P##_F_from_site(F_OFFSET(mom), EVENANDODD); \
  if(*rawmom == NULL)terminate(1); \
  *mom = QOP_##P##3_create_F_from_raw((MILCFLOAT **)(*rawmom),QOP_EVENODD); \
}

#define make_load_links_and_mom_field(P, MILCTYPE, MILCFLOAT) \
void \
load_##P##_links_and_mom_field(QOP_##P##3_GaugeField **links, \
   QOP_##P##3_Force **mom, MILCTYPE ***rawlinks, MILCTYPE ***rawmom, \
   su3_matrix *srclink, anti_hermitmat *srcmom) \
{ \
  *rawlinks = create_raw4_##P##_G_from_field(srclink, EVENANDODD); \
  if(*rawlinks == NULL)terminate(1); \
  *links = QOP_##P##3_create_G_from_raw((MILCFLOAT **)(*rawlinks),QOP_EVENODD); \
  *rawmom = create_raw4_##P##_F_from_field(srcmom, EVENANDODD); \
  if(*rawmom == NULL)terminate(1); \
  *mom = QOP_##P##3_create_F_from_raw((MILCFLOAT **)(*rawmom),QOP_EVENODD); \
}

#define make_load_from_site(P, T, TYPE, MILCTYPE, MILCFLOAT) \
void \
load_##P##_##T##_from_site( TYPE** qop, field_offset src, int parity) \
{ \
  MILCTYPE *raw; \
  raw = create_raw_##P##_##T##_from_site(src, parity); \
  if(raw == NULL)terminate(1); \
  *qop = QOP_##P##3_create_##T##_from_raw((MILCFLOAT *)raw, \
					  milc2qop_parity(parity)); \
  destroy_raw_##P##_##T(raw); raw = NULL; \
  return; \
}

/* Map color vector field from MILC field to QOP field */

#define make_load_from_field(P, T, TYPE, MILCTYPE, MILC_SRC_TYPE, MILCFLOAT) \
void \
load_##P##_##T##_from_field( TYPE** qop, MILC_SRC_TYPE *src, int parity) \
{ \
  MILCTYPE *raw; \
  raw = create_raw_##P##_##T##_from_field(src, parity); \
  if(raw == NULL)terminate(1); \
  *qop = QOP_##P##3_create_##T##_from_raw((MILCFLOAT *)raw, \
                                        milc2qop_parity(parity)); \
  destroy_raw_##P##_##T(raw); raw = NULL; \
  return; \
}

#define make_unload_links_and_mom_site(P, MILCTYPE, MILCFLOAT) \
void \
unload_##P##_links_and_mom_site(QOP_##P##3_GaugeField **links,  \
   QOP_##P##3_Force **mom, MILCTYPE ***rawlinks, MILCTYPE ***rawmom) \
{ \
  destroy_raw4_##P##_G (*rawlinks);   *rawlinks = NULL; \
  QOP_##P##3_destroy_G (*links);      *links = NULL; \
  QOP_##P##3_extract_F_to_raw((MILCFLOAT **)(*rawmom), *mom, QOP_EVENODD); \
  unload_raw4_##P##_F_to_site(F_OFFSET(mom), *rawmom, EVENANDODD); \
  destroy_raw4_##P##_F (*rawmom);   *rawmom = NULL; \
  QOP_##P##3_destroy_F (*mom);      *mom = NULL; \
}

#define make_unload_links_and_mom_field(P, MILCTYPE, MILCFLOAT) \
void \
unload_##P##_links_and_mom_field(su3_matrix *dstlink, anti_hermitmat *dstmom, \
   QOP_##P##3_GaugeField **links, QOP_##P##3_Force **mom, \
   MILCTYPE ***rawlinks, MILCTYPE ***rawmom) \
{ \
  destroy_raw4_##P##_G (*rawlinks);   *rawlinks = NULL; \
  QOP_##P##3_destroy_G (*links);      *links = NULL; \
  QOP_##P##3_extract_F_to_raw((MILCFLOAT **)(*rawmom), *mom, QOP_EVENODD); \
  unload_raw4_##P##_F_to_field(dstmom, *rawmom, EVENANDODD); \
  destroy_raw4_##P##_F (*rawmom);   *rawmom = NULL; \
  QOP_##P##3_destroy_F (*mom);      *mom = NULL; \
}

/* Map color vector from QOP field to site */

#define make_unload_to_site(P, T, TYPE, MILCTYPE, MILCFLOAT) \
void \
unload_##P##_##T##_to_site( field_offset dest, TYPE *qop, int parity){ \
  MILCTYPE *raw; \
  raw = create_raw_##P##_##T(); \
  QOP_##P##3_extract_##T##_to_raw((MILCFLOAT *)raw, qop, milc2qop_parity(parity)); \
  unload_raw_##P##_##T##_to_site(dest, raw, parity); \
  destroy_raw_##P##_##T(raw); raw = NULL; \
}

/* Map color vector from QOP field to MILC field */

#define make_unload_to_field(P, T, TYPE, MILCTYPE, MILC_DSTTYPE, MILCFLOAT) \
void \
unload_##P##_##T##_to_field( MILC_DSTTYPE *dest, TYPE *qop, int parity){ \
  MILCTYPE *raw; \
  raw = create_raw_##P##_##T(); \
  QOP_##P##3_extract_##T##_to_raw((MILCFLOAT *)raw, qop, milc2qop_parity(parity)); \
  unload_raw_##P##_##T##_to_field(dest, raw, parity); \
  destroy_raw_##P##_##T(raw); raw = NULL; \
}

/* Storage for raw gauge field */

make_create_raw4(F, G, fsu3_matrix);
make_create_raw4(D, G, dsu3_matrix);

make_destroy_raw4(F, G, fsu3_matrix);
make_destroy_raw4(D, G, dsu3_matrix);

/* Storage for raw gauge momentum */

make_create_raw4(F, F, fsu3_matrix);
make_create_raw4(D, F, dsu3_matrix);

make_destroy_raw4(F, F, fsu3_matrix);
make_destroy_raw4(D, F, dsu3_matrix);

/* Storage for raw su3 vector field */

make_create_raw(F, V, fsu3_vector);
make_create_raw(D, V, dsu3_vector);

make_destroy_raw(F, V, fsu3_vector);
make_destroy_raw(D, V, dsu3_vector);

/* Map gauge field from site to raw */

make_create_raw4_from_site(F, G, fsu3_matrix, su3_matrix);
make_create_raw4_from_site(D, G, dsu3_matrix, su3_matrix);

/* Map gauge field from field to raw */

make_create_raw4_from_field(F, G, fsu3_matrix, su3_matrix);
make_create_raw4_from_field(D, G, dsu3_matrix, su3_matrix);

/* Map gauge momentum from site to raw */

make_create_raw4_from_site(F, F, fsu3_matrix, anti_hermitmat);
make_create_raw4_from_site(D, F, dsu3_matrix, anti_hermitmat);

/* Map gauge momentum from field to raw */

make_create_raw4_from_field(F, F, fsu3_matrix, anti_hermitmat);
make_create_raw4_from_field(D, F, dsu3_matrix, anti_hermitmat);

/* Map color vector from site to raw */

make_create_raw_from_site(F, V, fsu3_vector, su3_vector);
make_create_raw_from_site(D, V, dsu3_vector, su3_vector);

/* Map color vector from field to raw */

make_create_raw_from_field(F, V, fsu3_vector, su3_vector);
make_create_raw_from_field(D, V, dsu3_vector, su3_vector);

/* Map gauge field from raw to site */

make_unload_raw4_to_site(F, G, su3_matrix, fsu3_matrix);
make_unload_raw4_to_site(D, G, su3_matrix, dsu3_matrix);

/* Map gauge field from raw to field */

make_unload_raw4_to_field(F, G, su3_matrix, fsu3_matrix);
make_unload_raw4_to_field(D, G, su3_matrix, dsu3_matrix);

/* Map gauge momentum from raw to site */

make_unload_raw4_to_site(F, F, anti_hermitmat, fsu3_matrix);
make_unload_raw4_to_site(D, F, anti_hermitmat, dsu3_matrix);

/* Map gauge momentum from raw to field */

make_unload_raw4_to_field(F, F, anti_hermitmat, fsu3_matrix);
make_unload_raw4_to_field(D, F, anti_hermitmat, dsu3_matrix);

/* Map color vector from raw to site */

make_unload_raw_to_site(F, V, su3_vector, fsu3_vector);
make_unload_raw_to_site(D, V, su3_vector, dsu3_vector);

/* Map color vector from raw to field */

make_unload_raw_to_field(F, V, su3_vector, fsu3_vector);
make_unload_raw_to_field(D, V, su3_vector, dsu3_vector);

/* Composite mapping */

make_load_links_and_mom_site(F, fsu3_matrix, float);
make_load_links_and_mom_site(D, dsu3_matrix, double);

make_load_links_and_mom_field(F, fsu3_matrix, float);
make_load_links_and_mom_field(D, dsu3_matrix, double);

make_unload_links_and_mom_site(F, fsu3_matrix, float);
make_unload_links_and_mom_site(D, dsu3_matrix, double);

make_unload_links_and_mom_field(F, fsu3_matrix, float);
make_unload_links_and_mom_field(D, dsu3_matrix, double);

make_load_from_site(F, V, QOP_F3_ColorVector, fsu3_vector, float);
make_load_from_site(D, V, QOP_D3_ColorVector, dsu3_vector, double);

make_load_from_field(F, V, QOP_F3_ColorVector, fsu3_vector, su3_vector ,float);
make_load_from_field(D, V, QOP_D3_ColorVector, dsu3_vector, su3_vector ,double);

make_unload_to_site(F, V, QOP_F3_ColorVector, fsu3_vector, float);
make_unload_to_site(D, V, QOP_D3_ColorVector, dsu3_vector, double);

make_unload_to_field(F, V, QOP_F3_ColorVector, fsu3_vector, su3_vector, float);
make_unload_to_field(D, V, QOP_D3_ColorVector, dsu3_vector, su3_vector, double);

/* map_milc_to_qop.c */

