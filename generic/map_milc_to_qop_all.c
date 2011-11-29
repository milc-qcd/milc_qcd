/*********************** map_milc_to_qop_all.c **************************/
/* Functions for mapping MILC data layouts to raw QOP layouts
   C. DeTar 4/29/07    add DiracFermion = wilson_fermion
   C. DeTar 5/15/07    split QDP from MILC/QCDOC
*/

   /* NOTE: This is an include file for map_milc_to_qopqdp.c and
      map_milc_to_qopmilc.c. */

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

#define copy_milc_to_F_D(d,s) p2f_wvec(d,s);
#define copy_milc_to_D_D(d,s) p2d_wvec(d,s);

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

#define copy_F_D_to_milc(d,s) f2p_wvec(d,s);
#define copy_D_D_to_milc(d,s) d2p_wvec(d,s);

void site_coords(int coords[4],site *s){
  coords[0] = s->x;
  coords[1] = s->y;
  coords[2] = s->z;
  coords[3] = s->t;
}

/* Map MILC site links to raw order */

#define make_create_raw4_from_site(P, T, RAWTYPE, MILC_SRCTYPE) \
RAWTYPE ** \
create_raw4_##P##_##T##_from_site(field_offset src, int milc_parity){ \
  int coords[4]; \
  int i,j,dir; \
  site *s; \
  RAWTYPE **raw; \
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

/* Map MILC field links to raw order */

#define make_create_raw4_from_field(P, T, RAWTYPE, MILC_SRCTYPE) \
RAWTYPE ** \
create_raw4_##P##_##T##_from_field(MILC_SRCTYPE *src, int milc_parity){ \
  int coords[4]; \
  int i,j,dir; \
  site *s; \
  RAWTYPE **raw = NULL; \
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

#define make_create_raw_from_site(P, T, RAWTYPE, MILC_SRCTYPE) \
RAWTYPE * \
create_raw_##P##_##T##_from_site(field_offset src, int milc_parity){ \
  int coords[4]; \
  int i,j; \
  site *s; \
  RAWTYPE *raw; \
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
    tmp = (MILC_SRCTYPE *)F_PT(s, src); \
    copy_milc_to_##P##_##T(raw + j, tmp); \
  } \
  return raw; \
}

/* Map MILC field to raw */

#define make_create_raw_from_field(P, T, RAWTYPE, MILC_SRCTYPE) \
RAWTYPE * \
create_raw_##P##_##T##_from_field(MILC_SRCTYPE *src, int milc_parity){ \
  int coords[4]; \
  int i,j; \
  site *s; \
  RAWTYPE *raw; \
  raw = create_raw_##P##_##T(); \
  if(raw == NULL)return NULL; \
  FORSOMEPARITY(i,s,milc_parity){ \
    site_coords(coords,s); \
    if(QOP_node_number_raw(coords) != this_node){ \
      printf("create_raw_from_field: incompatible layout\n"); \
      return NULL; \
    } \
    j = QOP_node_index_raw_##T (coords, milc2qop_parity(milc_parity)); \
    copy_milc_to_##P##_##T(raw + j, src + i); \
  } \
  return raw; \
}

/* Map MILC gauge field in site structure to QOP through raw type */

#define make_create_from_site4(P, T, QOPTYPE, RAWTYPE, MILCFLOAT) \
QOPTYPE * \
create_##P##_##T##_from_site4(field_offset src, int milc_parity){ \
  RAWTYPE **raw; \
  QOPTYPE *qop; \
  raw = create_raw4_##P##_##T##_from_site(src, milc_parity); \
  if(raw == NULL)terminate(1); \
  qop = QOP_##P##3_create_##T##_from_raw((MILCFLOAT **)raw, \
           milc2qop_parity(milc_parity)); \
  destroy_raw4_##P##_##T(raw); raw = NULL; \
  return qop; \
}

/* Map MILC field to QOP through raw type */

#define make_create_from_field4(P, T, QOPTYPE, RAWTYPE, MILC_SRCTYPE, MILCFLOAT) \
QOPTYPE * \
create_##P##_##T##_from_field4(MILC_SRCTYPE *src, int milc_parity){ \
  RAWTYPE **raw; \
  QOPTYPE *qop; \
  raw = create_raw4_##P##_##T##_from_field(src, milc_parity); \
  if(raw == NULL)terminate(1); \
  qop = QOP_##P##3_create_##T##_from_raw((MILCFLOAT **)raw, \
           milc2qop_parity(milc_parity)); \
  destroy_raw4_##P##_##T(raw); raw = NULL; \
  return qop; \
}

/* Map MILC site to QOP through raw type */

#define make_create_from_site(P, T, QOPTYPE, RAWTYPE, MILCFLOAT) \
QOPTYPE * \
create_##P##_##T##_from_site(field_offset src, int milc_parity){ \
  RAWTYPE *raw; \
  QOPTYPE *qop; \
  raw = create_raw_##P##_##T##_from_site(src, milc_parity); \
  if(raw == NULL)terminate(1); \
  qop = QOP_##P##3_create_##T##_from_raw((MILCFLOAT *)raw, \
           milc2qop_parity(milc_parity)); \
  destroy_raw_##P##_##T(raw); raw = NULL; \
  return qop; \
}

/* Map MILC field to QOP through raw type */

#define make_create_from_field(P, T, QOPTYPE, RAWTYPE, MILC_SRCTYPE, MILCFLOAT) \
QOPTYPE * \
create_##P##_##T##_from_field(MILC_SRCTYPE *src, int milc_parity){ \
  RAWTYPE *raw; \
  QOPTYPE *qop; \
  raw = create_raw_##P##_##T##_from_field(src, milc_parity); \
  if(raw == NULL)terminate(1); \
  qop = QOP_##P##3_create_##T##_from_raw((MILCFLOAT *)raw, \
           milc2qop_parity(milc_parity)); \
  destroy_raw_##P##_##T(raw); raw = NULL; \
  return qop; \
}

/* Map raw links to MILC site structure */

#define make_unload_raw4_to_site(P, T, MILC_DSTTYPE, RAWTYPE) \
void \
unload_raw4_##P##_##T##_to_site(field_offset dest, RAWTYPE *raw[], \
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

#define make_unload_raw4_to_field(P, T, MILC_DSTTYPE, RAWTYPE) \
void \
unload_raw4_##P##_##T##_to_field(MILC_DSTTYPE *dest, RAWTYPE *raw[], \
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

#define make_unload_raw_to_site(P, T, MILC_DSTTYPE, RAWTYPE) \
void \
unload_raw_##P##_##T##_to_site(field_offset dest, RAWTYPE *raw, \
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

#define make_unload_raw_to_field(P, T, MILC_DSTTYPE, RAWTYPE) \
void \
unload_raw_##P##_##T##_to_field(MILC_DSTTYPE *dest, RAWTYPE *raw, \
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

/* Map MILC gauge field in site structure to QOP through raw type */

#define make_unload_to_site4(P, T, QOPTYPE, RAWTYPE, MILCFLOAT) \
void \
unload_##P##_##T##_to_site4(field_offset dest, QOPTYPE *qop, int milc_parity){ \
  RAWTYPE **raw; \
  raw = create_raw4_##P##_##T(); \
  if(raw == NULL)terminate(1); \
  QOP_##P##3_extract_##T##_to_raw((MILCFLOAT **)raw, qop, \
           milc2qop_parity(milc_parity)); \
  unload_raw4_##P##_##T##_to_site(dest, raw, milc_parity); \
  destroy_raw4_##P##_##T(raw); raw = NULL; \
  return; \
}

/* Map MILC field to QOP through raw type */

#define make_unload_to_field4(P, T, QOPTYPE, RAWTYPE, MILC_DSTTYPE, MILCFLOAT) \
void \
 unload_##P##_##T##_to_field4(MILC_DSTTYPE *dest, QOPTYPE *qop, int milc_parity){ \
  RAWTYPE **raw; \
  raw = create_raw4_##P##_##T(); \
  if(raw == NULL)terminate(1); \
  QOP_##P##3_extract_##T##_to_raw((MILCFLOAT **)raw, qop, \
           milc2qop_parity(milc_parity)); \
  unload_raw4_##P##_##T##_to_field(dest, raw, milc_parity); \
  destroy_raw4_##P##_##T(raw); raw = NULL; \
  return; \
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

#if 0
/* Generate QOP asqtad links from gauge field in site structure */

#define make_create_L_from_site_gauge(P, TYPE, MILCTYPE, MILCFLOAT) \
TYPE * \
create_##P##_L_from_site_gauge( QOP_info_t *info, \
    QOP_asqtad_coeffs_t *coeffs, field_offset src, int parity) \
{ \
  MILCTYPE **raw; \
  TYPE *qop; \
  QOP_##P##3_GaugeField *gauge; \
  raw = create_raw4_##P##_G_from_site(src, parity); \
  if(raw == NULL)terminate(1); \
  gauge = QOP_##P##3_create_G_from_raw((MILCFLOAT **)raw, \
				  milc2qop_parity(parity)); \
  destroy_raw4_##P##_G(raw); raw = NULL; \
  qop = QOP_##P##3_asqtad_create_L_from_G(info, coeffs, gauge); \
  QOP_##P##3_destroy_G(gauge); gauge = NULL; \
  return qop; \
}
#endif

/* Map preconstructed fat and long gauge fields from MILC site to QOP field */

#define make_create_L_from_sites(P, TYPE, MILCTYPE, MILCFLOAT) \
TYPE * \
create_##P##_L_from_sites( field_offset fat, field_offset lng, \
   int parity) \
{ \
  MILCTYPE **rawfat; \
  MILCTYPE **rawlng; \
  TYPE *qop; \
  rawfat = create_raw4_##P##_G_from_site(fat, parity); \
  if(rawfat == NULL)terminate(1); \
  rawlng = create_raw4_##P##_G_from_site(lng, parity); \
  if(rawlng == NULL)terminate(1); \
  qop = QOP_##P##3_asqtad_create_L_from_raw((MILCFLOAT **)rawfat, \
          (MILCFLOAT **)rawlng, milc2qop_parity(parity)); \
  destroy_raw4_##P##_G(rawfat); rawfat = NULL; \
  destroy_raw4_##P##_G(rawlng); rawlng = NULL; \
  return qop; \
}

/* Map preconstructed fat and long gauge fields from MILC field to QOP field */

#define make_create_L_from_fields(P, TYPE, MILCTYPE, MILC_SRC_TYPE, MILCFLOAT) \
TYPE * \
create_##P##_L_from_fields( MILC_SRC_TYPE *fat, MILC_SRC_TYPE *lng, \
   int parity) \
{ \
  MILCTYPE **rawfat; \
  MILCTYPE **rawlng; \
  TYPE *qop; \
  rawfat = create_raw4_##P##_G_from_field(fat, parity); \
  if(rawfat == NULL)terminate(1); \
  rawlng = create_raw4_##P##_G_from_field(lng, parity); \
  if(rawlng == NULL)terminate(1); \
  qop = QOP_##P##3_asqtad_create_L_from_raw((MILCFLOAT **)rawfat, \
           (MILCFLOAT **)rawlng, milc2qop_parity(parity)); \
  destroy_raw4_##P##_G(rawfat); rawfat = NULL; \
  destroy_raw4_##P##_G(rawlng); rawlng = NULL; \
  return qop; \
}

/* Extract MILC fat link and long link fields from QOP asqtad links */

#define make_unload_L_to_fields(P, TYPE, MILCTYPE, MILC_DST_TYPE, MILCFLOAT) \
void \
unload_##P##_L_to_fields( MILC_DST_TYPE *fat, MILC_DST_TYPE *lng, TYPE* qop, \
   int parity) \
{ \
  MILCTYPE **rawfat; \
  MILCTYPE **rawlng; \
  rawfat  = create_raw4_##P##_G(); \
  if(rawfat == NULL)terminate(1); \
  rawlng = create_raw4_##P##_G(); \
  if(rawlng == NULL)terminate(1); \
  QOP_##P##3_asqtad_extract_L_to_raw((MILCFLOAT **)rawfat, \
          (MILCFLOAT **)rawlng, qop, milc2qop_parity(parity)); \
  unload_raw4_##P##_G_to_field(fat, rawfat, parity); \
  if(lng != NULL) \
  unload_raw4_##P##_G_to_field(lng, rawlng, parity); \
  destroy_raw4_##P##_G(rawfat); rawfat = NULL; \
  destroy_raw4_##P##_G(rawlng); rawlng = NULL; \
  return; \
}


/* Extract MILC fat link and long link fields from QOP HISQ links */

#define make_unload_hisq_L_to_fields(P, TYPE, MILCTYPE, MILC_DST_TYPE, MILCFLOAT) \
void \
unload_##P##_hisq_L_to_fields( MILC_DST_TYPE *fat, MILC_DST_TYPE *lng, TYPE* qop, \
   int parity) \
{ \
  MILCTYPE **rawfat; \
  MILCTYPE **rawlng; \
  rawfat  = create_raw4_##P##_G(); \
  if(rawfat == NULL)terminate(1); \
  rawlng = create_raw4_##P##_G(); \
  if(rawlng == NULL)terminate(1); \
  QOP_##P##3_hisq_extract_L_to_raw((MILCFLOAT **)rawfat, \
          (MILCFLOAT **)rawlng, qop, milc2qop_parity(parity)); \
  unload_raw4_##P##_G_to_field(fat, rawfat, parity); \
  if(lng != NULL) \
  unload_raw4_##P##_G_to_field(lng, rawlng, parity); \
  destroy_raw4_##P##_G(rawfat); rawfat = NULL; \
  destroy_raw4_##P##_G(rawlng); rawlng = NULL; \
  return; \
}


/* Map MILC clover term to raw QOP clover */


#define DI(x) (-(1 - (x))/2.)
#define TR(x) ((x)/2.)

#define make_map_milc_clov_to_qop_raw(P, MILCFLOAT) \
void map_milc_clov_to_qop_raw_##P(MILCFLOAT *raw_clov, clover *milc_clov ){\
  int i,j,c;\
  int coords[4]; \
  site *s;\
  MILCFLOAT *r;\
\
  FORALLSITES(i,s){\
\
    site_coords(coords,s);\
    if(QOP_node_number_raw(coords) != this_node){\
      printf("map_milc_clov_to_qop_raw: QOP layout is incompatible with MILC layout\n");\
      terminate(1);\
    }\
    j = QOP_node_index_raw_D(coords, QOP_EVENODD);\
    r = raw_clov + 72*j;\
\
    for(c = 0; c < 3; c++){\
      r[2*c]   = DI(milc_clov->clov_diag[i].di[0][c]);      /* c0 c0 */\
      r[2*c+1] = DI(milc_clov->clov_diag[i].di[0][c+3]);    /* c1 c1 */\
    }\
      \
    r += 6;\
\
    r[ 0] = TR( milc_clov->clov[i].tr[0][ 3].real);         /* 01 00 */\
    r[ 1] = TR( milc_clov->clov[i].tr[0][ 3].imag);         /* 01 00 */\
    r[ 2] = TR( milc_clov->clov[i].tr[0][ 0].real);         /* 10 00 */\
    r[ 3] = TR( milc_clov->clov[i].tr[0][ 0].imag);         /* 10 00 */\
    r[ 4] = TR( milc_clov->clov[i].tr[0][ 6].real);         /* 11 00 */\
    r[ 5] = TR( milc_clov->clov[i].tr[0][ 6].imag);         /* 11 00 */\
    r[ 6] = TR( milc_clov->clov[i].tr[0][ 1].real);         /* 20 00 */\
    r[ 7] = TR( milc_clov->clov[i].tr[0][ 1].imag);         /* 20 00 */\
    r[ 8] = TR( milc_clov->clov[i].tr[0][10].real);         /* 21 00 */\
    r[ 9] = TR( milc_clov->clov[i].tr[0][10].imag);         /* 21 00 */\
\
    r[10] = TR( milc_clov->clov[i].tr[0][ 4].real);         /* 10 01 */\
    r[11] = TR(-milc_clov->clov[i].tr[0][ 4].imag);         /* 10 01 */\
    r[12] = TR( milc_clov->clov[i].tr[0][ 9].real);         /* 11 01 */\
    r[13] = TR( milc_clov->clov[i].tr[0][ 9].imag);         /* 11 01 */\
    r[14] = TR( milc_clov->clov[i].tr[0][ 5].real);         /* 20 01 */\
    r[15] = TR(-milc_clov->clov[i].tr[0][ 5].imag);         /* 20 01 */\
    r[16] = TR( milc_clov->clov[i].tr[0][13].real);         /* 21 01 */\
    r[17] = TR( milc_clov->clov[i].tr[0][13].imag);         /* 21 01 */\
\
    r[18] = TR( milc_clov->clov[i].tr[0][ 7].real);         /* 11 10 */\
    r[19] = TR( milc_clov->clov[i].tr[0][ 7].imag);         /* 11 10 */\
    r[20] = TR( milc_clov->clov[i].tr[0][ 2].real);         /* 20 10 */\
    r[21] = TR( milc_clov->clov[i].tr[0][ 2].imag);         /* 20 10 */\
    r[22] = TR( milc_clov->clov[i].tr[0][11].real);         /* 21 10 */\
    r[23] = TR( milc_clov->clov[i].tr[0][11].imag);         /* 21 10 */\
\
    r[24] = TR( milc_clov->clov[i].tr[0][ 8].real);         /* 20 11 */\
    r[25] = TR(-milc_clov->clov[i].tr[0][ 8].imag);         /* 20 11 */\
    r[26] = TR( milc_clov->clov[i].tr[0][14].real);         /* 21 11 */\
    r[27] = TR( milc_clov->clov[i].tr[0][14].imag);         /* 21 11 */\
\
    r[28] = TR( milc_clov->clov[i].tr[0][12].real);         /* 21 20 */\
    r[29] = TR( milc_clov->clov[i].tr[0][12].imag);         /* 21 20 */\
\
    r += 30;\
\
    for(c = 0; c < 3; c++){\
      r[2*c]   = DI(milc_clov->clov_diag[i].di[1][c]);      /* c2 c2 */\
      r[2*c+1] = DI(milc_clov->clov_diag[i].di[1][c+3]);    /* c3 c3 */\
    }\
\
    r += 6;\
\
    r[ 0] = TR( milc_clov->clov[i].tr[1][ 3].real);         /* 03 02 */\
    r[ 1] = TR( milc_clov->clov[i].tr[1][ 3].imag);         /* 03 02 */\
    r[ 2] = TR( milc_clov->clov[i].tr[1][ 0].real);         /* 12 02 */\
    r[ 3] = TR( milc_clov->clov[i].tr[1][ 0].imag);         /* 12 02 */\
    r[ 4] = TR( milc_clov->clov[i].tr[1][ 6].real);         /* 13 02 */\
    r[ 5] = TR( milc_clov->clov[i].tr[1][ 6].imag);         /* 13 02 */\
    r[ 6] = TR( milc_clov->clov[i].tr[1][ 1].real);         /* 22 02 */\
    r[ 7] = TR( milc_clov->clov[i].tr[1][ 1].imag);         /* 22 02 */\
    r[ 8] = TR( milc_clov->clov[i].tr[1][10].real);         /* 23 02 */\
    r[ 9] = TR( milc_clov->clov[i].tr[1][10].imag);         /* 23 02 */\
\
    r[10] = TR( milc_clov->clov[i].tr[1][ 4].real);         /* 12 03 */\
    r[11] = TR(-milc_clov->clov[i].tr[1][ 4].imag);         /* 12 03 */\
    r[12] = TR( milc_clov->clov[i].tr[1][ 9].real);         /* 13 03 */\
    r[13] = TR( milc_clov->clov[i].tr[1][ 9].imag);         /* 13 03 */\
    r[14] = TR( milc_clov->clov[i].tr[1][ 5].real);         /* 22 03 */\
    r[15] = TR(-milc_clov->clov[i].tr[1][ 5].imag);         /* 22 03 */\
    r[16] = TR( milc_clov->clov[i].tr[1][13].real);         /* 23 03 */\
    r[17] = TR( milc_clov->clov[i].tr[1][13].imag);         /* 23 03 */\
\
    r[18] = TR( milc_clov->clov[i].tr[1][ 7].real);         /* 13 12 */\
    r[19] = TR( milc_clov->clov[i].tr[1][ 7].imag);         /* 13 12 */\
    r[20] = TR( milc_clov->clov[i].tr[1][ 2].real);         /* 22 12 */\
    r[21] = TR( milc_clov->clov[i].tr[1][ 2].imag);         /* 22 12 */\
    r[22] = TR( milc_clov->clov[i].tr[1][11].real);         /* 23 12 */\
    r[23] = TR( milc_clov->clov[i].tr[1][11].imag);         /* 23 12 */\
\
    r[24] = TR( milc_clov->clov[i].tr[1][ 8].real);         /* 22 13 */\
    r[25] = TR(-milc_clov->clov[i].tr[1][ 8].imag);         /* 22 13 */\
    r[26] = TR( milc_clov->clov[i].tr[1][14].real);         /* 23 13 */\
    r[27] = TR( milc_clov->clov[i].tr[1][14].imag);         /* 23 13 */\
\
    r[28] = TR( milc_clov->clov[i].tr[1][12].real);         /* 23 22 */\
    r[29] = TR( milc_clov->clov[i].tr[1][12].imag);         /* 23 22 */\
\
    r += 30;\
  }\
}

/* map_milc_to_qop_all.c */

