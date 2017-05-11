/*********************** map_milc_to_qphixj_all.c **************************/
/* Functions for mapping MILC data layouts to raw QPHIXJ layouts
   C. DeTar 12/23/15    add DiracFermion = wilson_fermion
*/

   /* NOTE: This is an include file for map_milc_to_qphixjqdp.c and
      map_milc_to_qphixjmilc.c. */

/* Temporary redefines. These definitions belong in the Qphixj interface. */

static int QPHIXJ_node_number_raw(int coords[]){
  return node_number(coords[0],coords[1],coords[2],coords[3]);
}

static int QPHIXJ_node_index_raw_G(int coords[], int milc_parity){
  return node_index(coords[0],coords[1],coords[2],coords[3]);
}

static int QPHIXJ_node_index_raw_F(int coords[], int milc_parity){
  return node_index(coords[0],coords[1],coords[2],coords[3]);
}

static int QPHIXJ_node_index_raw_V(int coords[], int milc_parity){
  return node_index(coords[0],coords[1],coords[2],coords[3]);
}

static int QPHIXJ_node_index_raw_D(int coords[], int milc_parity){
  return node_index(coords[0],coords[1],coords[2],coords[3]);
}

/* Create empty raw links */
/* Note, this version assumes four contiguous matrices per site */

#define make_create_raw4(P, T, MILCTYPE) \
MILCTYPE * \
create_qphixj_raw4_##P##_##T (void){ \
 MILCTYPE *raw = NULL; \
 int i; \
 raw = (MILCTYPE *)malloc(4*sites_on_node*sizeof(MILCTYPE)); \
 if(raw == NULL){                              \
   printf("create4_qphixj_raw: No room\n");     \
   return NULL;                                \
 } \
 FORALLFIELDSITES_OMP(i, ){ \
 memset(raw+4*i, 0, 4*sizeof(MILCTYPE));     \
 } END_LOOP_OMP; \
 return raw; \
}

/* Destroy raw links */
/* Note, this version assumes four contiguous matrices per site */

#define make_destroy_raw4(P, T, MILCTYPE) \
void \
destroy_qphixj_raw4_##P##_##T (MILCTYPE *raw){ \
  free(raw); \
}

/* Create empty raw field */

#define make_create_raw(P, T, MILCTYPE) \
MILCTYPE * \
create_qphixj_raw_##P##_##T(void){ \
  MILCTYPE *raw = NULL; \
  int i; \
  raw = (MILCTYPE *)malloc(sites_on_node*sizeof(MILCTYPE)); \
  if(raw == NULL){ \
    printf("create_qphixj_raw: No room\n"); \
    return NULL; \
  } \
  FORALLFIELDSITES_OMP(i, ){ \
    memset(raw+i, 0, sizeof(MILCTYPE));     \
  } END_LOOP_OMP; \
  return raw; \
}

/* Destroy raw field */

#define make_destroy_raw(P, T, MILCTYPE) \
void \
destroy_qphixj_raw_##P##_##T (MILCTYPE *raw){ \
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

static void 
site_coords(int coords[4],site *s){
  coords[0] = s->x;
  coords[1] = s->y;
  coords[2] = s->z;
  coords[3] = s->t;
}

static int layout_compatibility_checked = 0;

static void check_layout_compatibility(void){
  int i;
  site *s;
  int coords[4];

  if(layout_compatibility_checked)return;

  FORALLSITES_OMP(i,s,private(coords)){
    site_coords(coords,s);
    if(QPHIXJ_node_number_raw(coords) != this_node){
      printf("ERROR on node(%d): incompatible QPHIXJ_node_number in layout\n", this_node);
      terminate(1);
    }
  } END_LOOP_OMP;
  layout_compatibility_checked = 1;
  return;
}

/* Note, there is actually no dependence on parity here */

#define make_create_hash(T)				\
static int hash_created_##T = 0;		\
static int *hash_m2q_##T;						\
static int *hash_q2m_##T;						\
								\
static void							\
 create_hash_##T(int milc_parity){				\
  site *s;						\
  int i, j, coords[4];					\
  if(!hash_created_##T){					\
     hash_m2q_##T = (int *)malloc(sizeof(int)*sites_on_node); \
     hash_q2m_##T = (int *)malloc(sizeof(int)*sites_on_node); \
     FORALLSITES_OMP(i,s,private(j,coords)){					\
       site_coords(coords,s);						\
       j = QPHIXJ_node_index_raw_##T(coords, milc2qphixj_parity(milc_parity)); \
       hash_m2q_##T[i] = j;							\
       hash_q2m_##T[j] = i;							\
     } END_LOOP_OMP;							\
     hash_created_##T = 1;						\
  }									\
}									\

/* Map MILC site links to raw order */
/* Note that unlike QOP we want four contiguous SU(3) matrices per site */

#define make_create_raw4_from_site(P, T, RAWTYPE, MILC_SRCTYPE) \
RAWTYPE * \
create_qphixj_raw4_##P##_##T##_from_site(field_offset src, int milc_parity){ \
  int i,j,dir,status=0; \
  site *s; \
  RAWTYPE *raw; \
  MILC_SRCTYPE *tmp; \
  raw = create_qphixj_raw4_##P##_##T (); \
  if(raw == NULL)return NULL; \
  check_layout_compatibility(); \
  create_hash_##T(milc_parity);					\
  FORSOMEPARITY_OMP(i,s,milc_parity,private(j,dir,tmp)){	\
    j = hash_m2q_##T[i];	   \
    if(i == 0)node0_printf("Hash for %d is %d\n",i,j);  \
    FORALLUPDIR(dir){ \
      tmp = (MILC_SRCTYPE *)F_PT(s, src); \
      copy_milc_to_##P##_##T(raw + 4*j + dir, tmp + dir); \
    } \
  } END_LOOP_OMP; \
  if(status == 1){ \
    destroy_qphixj_raw4_##P##_##T (raw); \
    return NULL; \
  } else { \
    return raw;\
  } \
}

/* Map MILC field links to raw order */
/* Note that unlike QOP we want four contiguous SU(3) matrices per site */

#define make_create_raw4_from_field(P, T, RAWTYPE, MILC_SRCTYPE) \
RAWTYPE * \
create_qphixj_raw4_##P##_##T##_from_field(MILC_SRCTYPE *src, int milc_parity){ \
  int i,j,dir,status=0;	 \
  site *s; \
  RAWTYPE *raw = NULL; \
  MILC_SRCTYPE *tmp; \
  raw = create_qphixj_raw4_##P##_##T (); \
  if(raw == NULL)return NULL; \
  check_layout_compatibility(); \
  create_hash_##T(milc_parity);					\
  FORSOMEPARITY_OMP(i,s,milc_parity,private(j,dir,tmp)){	\
    j = hash_m2q_##T[i];	   \
    FORALLUPDIR(dir){ \
      tmp = src + 4*i; \
      copy_milc_to_##P##_##T(raw + 4*j + dir, tmp + dir); \
    } \
  } END_LOOP_OMP; \
  if(status == 1){ \
    destroy_qphixj_raw4_##P##_##T (raw); \
    return NULL; \
  } else { \
    return raw;\
  } \
}

/* Map MILC site field to raw */

#define make_create_raw_from_site(P, T, RAWTYPE, MILC_SRCTYPE) \
RAWTYPE * \
create_qphixj_raw_##P##_##T##_from_site(field_offset src, int milc_parity){ \
  int i,j,status=0; \
  site *s; \
  RAWTYPE *raw; \
  MILC_SRCTYPE *tmp; \
  raw = create_qphixj_raw_##P##_##T(); \
  if(raw == NULL)return NULL; \
  check_layout_compatibility(); \
  create_hash_##T(milc_parity);					\
  FORSOMEPARITY_OMP(i,s,milc_parity,private(j,tmp)){	\
    j = hash_m2q_##T[i];	   \
    tmp = (MILC_SRCTYPE *)F_PT(s, src); \
    copy_milc_to_##P##_##T(raw + j, tmp); \
  } END_LOOP_OMP; \
  if(status == 1){ \
    destroy_qphixj_raw_##P##_##T (raw); \
    return NULL; \
  } else { \
    return raw;\
  } \
}

/* Map MILC field to raw */

#define make_create_raw_from_field(P, T, RAWTYPE, MILC_SRCTYPE) \
RAWTYPE * \
create_qphixj_raw_##P##_##T##_from_field(MILC_SRCTYPE *src, int milc_parity){ \
  int i,j,status=0; \
  site *s; \
  RAWTYPE *raw; \
  raw = create_qphixj_raw_##P##_##T(); \
  if(raw == NULL)return NULL; \
  check_layout_compatibility(); \
  create_hash_##T(milc_parity);					\
  FORSOMEFIELDPARITY_OMP(i,milc_parity,private(j)){	\
    j = hash_m2q_##T[i];	   \
    copy_milc_to_##P##_##T(raw + j, src + i); \
  } END_LOOP_OMP; \
  if(status == 1){ \
    destroy_qphixj_raw_##P##_##T (raw); \
    return NULL; \
  } else { \
    return raw;\
  } \
}

/* Map MILC gauge field in site structure to QPHIXJ through raw type */
/* Note, this version assumes four contiguous matrices per site */

#define make_create_from_site4(P, T, QPHIXJTYPE, RAWTYPE, MILCFLOAT) \
QPHIXJTYPE * \
create_qphixj_##P##_##T##_from_site4(field_offset src, int milc_parity){ \
  RAWTYPE *raw; \
  QPHIXJTYPE *qphixj; \
  raw = create_qphixj_raw4_##P##_##T##_from_site(src, milc_parity); \
  if(raw == NULL)terminate(1); \
  qphixj = QPHIXJ_##P##3_create_##T##_from_raw((MILCFLOAT *)raw, \
           milc2qphixj_parity(milc_parity)); \
  destroy_qphixj_raw4_##P##_##T(raw); raw = NULL; \
  return qphixj; \
}

/* Map MILC field to QPHIXJ through raw type */
/* Note, this version assumes four contiguous matrices per site */

#define make_create_from_field4(P, T, QPHIXJTYPE, RAWTYPE, MILC_SRCTYPE, MILCFLOAT) \
QPHIXJTYPE * \
create_qphixj_##P##_##T##_from_field4(MILC_SRCTYPE *src, int milc_parity){ \
  RAWTYPE *raw; \
  QPHIXJTYPE *qphixj; \
  raw = create_qphixj_raw4_##P##_##T##_from_field(src, milc_parity); \
  if(raw == NULL)terminate(1); \
  qphixj = QPHIXJ_##P##3_create_##T##_from_raw((MILCFLOAT *)raw, \
           milc2qphixj_parity(milc_parity)); \
  destroy_qphixj_raw4_##P##_##T(raw); raw = NULL; \
  return qphixj; \
}

/* Map MILC site to QPHIXJ through raw type */

#define make_create_from_site(P, T, QPHIXJTYPE, RAWTYPE, MILCFLOAT) \
QPHIXJTYPE * \
create_qphixj_##P##_##T##_from_site(field_offset src, int milc_parity){ \
  RAWTYPE *raw; \
  QPHIXJTYPE *qphixj; \
  raw = create_qphixj_raw_##P##_##T##_from_site(src, milc_parity); \
  if(raw == NULL)terminate(1); \
  qphixj = QPHIXJ_##P##3_create_##T##_from_raw((MILCFLOAT *)raw, \
           milc2qphixj_parity(milc_parity)); \
  destroy_qphixj_raw_##P##_##T(raw); raw = NULL; \
  return qphixj; \
}

/* Map MILC field to QPHIXJ through raw type */

#define make_create_from_field(P, T, QPHIXJTYPE, RAWTYPE, MILC_SRCTYPE, MILCFLOAT) \
QPHIXJTYPE * \
create_qphixj_##P##_##T##_from_field(MILC_SRCTYPE *src, int milc_parity){ \
  RAWTYPE *raw; \
  QPHIXJTYPE *qphixj; \
  raw = create_qphixj_raw_##P##_##T##_from_field(src, milc_parity); \
  if(raw == NULL)terminate(1); \
  qphixj = QPHIXJ_##P##3_create_##T##_from_raw((MILCFLOAT *)raw, \
           milc2qphixj_parity(milc_parity)); \
  destroy_qphixj_raw_##P##_##T(raw); raw = NULL; \
  return qphixj; \
}

/* Map raw links to MILC site structure */
/* Note that unlike QOP we want four contiguous SU(3) matrices per site */

#define make_unload_raw4_to_site(P, T, MILC_DSTTYPE, RAWTYPE) \
void \
unload_qphixj_raw4_##P##_##T##_to_site(field_offset dest, RAWTYPE *raw, \
         int milc_parity){ \
  int i,j,dir,status=0;	 \
  site *s; \
  MILC_DSTTYPE *tmp; \
  check_layout_compatibility(); \
  create_hash_##T(milc_parity);					\
  FORSOMEPARITY_OMP(i,s,milc_parity,private(dir,j,tmp)){	\
    j = hash_q2m_##T[i];	   \
    FORALLUPDIR(dir){ \
      tmp = (MILC_DSTTYPE *)F_PT(s, dest); \
      copy_##P##_##T##_to_milc(&tmp[dir], raw + 4*j + dir); \
    } END_LOOP_OMP; \
    if(status==1)terminate(1); \
  } \
}

/* Note that unlike QOP we want four contiguous SU(3) matrices per site */

#define make_unload_raw4_to_field(P, T, MILC_DSTTYPE, RAWTYPE) \
void \
unload_qphixj_raw4_##P##_##T##_to_field(MILC_DSTTYPE *dest, RAWTYPE *raw, \
         int milc_parity){ \
  int i,j,dir,status=0;	 \
  site *s; \
  MILC_DSTTYPE *tmp; \
  check_layout_compatibility(); \
  create_hash_##T(milc_parity);					\
  FORSOMEPARITY_OMP(i,s,milc_parity,private(j,dir,tmp)){	\
    j = hash_q2m_##T[i];	   \
    FORALLUPDIR(dir){ \
      tmp = dest + 4*i; \
      copy_##P##_##T##_to_milc(&tmp[dir], raw + 4*j + dir); \
    } \
  } END_LOOP_OMP; \
  if(status==1)terminate(1); \
}

#define make_unload_raw_to_site(P, T, MILC_DSTTYPE, RAWTYPE) \
void \
unload_qphixj_raw_##P##_##T##_to_site(field_offset dest, RAWTYPE *raw, \
       int milc_parity){ \
  int i,j,status=0; \
  site *s; \
  check_layout_compatibility(); \
  create_hash_##T(milc_parity);					\
  FORSOMEPARITY_OMP(i,s,milc_parity,private(j)){	\
    j = hash_q2m_##T[i];	   \
    copy_##P##_##T##_to_milc((MILC_DSTTYPE *)F_PT(s,dest), raw + j); \
  } END_LOOP_OMP; \
  if(status==1)terminate(1); \
}

#define make_unload_raw_to_field(P, T, MILC_DSTTYPE, RAWTYPE) \
void \
unload_qphixj_raw_##P##_##T##_to_field(MILC_DSTTYPE *dest, RAWTYPE *raw, \
       int milc_parity){ \
  int i,j,status=0;  \
  site *s; \
  check_layout_compatibility(); \
  create_hash_##T(milc_parity);					\
  FORSOMEPARITY_OMP(i,s,milc_parity,private(j)){	\
    j = hash_q2m_##T[i];	   \
    copy_##P##_##T##_to_milc(dest + i, raw + j); \
  } END_LOOP_OMP; \
  if(status==1)terminate(1); \
}

/* Map MILC gauge field in site structure to QPHIXJ through raw type */

#define make_unload_to_site4(P, T, QPHIXJTYPE, RAWTYPE, MILCFLOAT) \
void \
unload_qphixj_##P##_##T##_to_site4(field_offset dest, QPHIXJTYPE *qphixj, int milc_parity){ \
  RAWTYPE **raw; \
  raw = create_qphixj_raw4_##P##_##T(); \
  if(raw == NULL)terminate(1); \
  QPHIXJ_##P##3_extract_##T##_to_raw((MILCFLOAT *)raw, qphixj, \
           milc2qphixj_parity(milc_parity)); \
  unload_qphixj_raw4_##P##_##T##_to_site(dest, raw, milc_parity); \
  destroy_qphixj_raw4_##P##_##T(raw); raw = NULL; \
  return; \
}

/* Map MILC field to QPHIXJ through raw type */

#define make_unload_to_field4(P, T, QPHIXJTYPE, RAWTYPE, MILC_DSTTYPE, MILCFLOAT) \
void \
 unload_qphixj_##P##_##T##_to_field4(MILC_DSTTYPE *dest, QPHIXJTYPE *qphixj, int milc_parity){ \
  RAWTYPE *raw; \
  raw = create_qphixj_raw4_##P##_##T(); \
  if(raw == NULL)terminate(1); \
  QPHIXJ_##P##3_extract_##T##_to_raw((MILCFLOAT *)raw, qphixj, \
           milc2qphixj_parity(milc_parity)); \
  unload_qphixj_raw4_##P##_##T##_to_field(dest, raw, milc_parity); \
  destroy_qphixj_raw4_##P##_##T(raw); raw = NULL; \
  return; \
}

/* Map color vector from QPHIXJ field to site */

#define make_unload_to_site(P, T, TYPE, MILCTYPE, MILCFLOAT) \
void \
unload_qphixj_##P##_##T##_to_site( field_offset dest, TYPE *qphixj, int parity){ \
  MILCTYPE *raw; \
  raw = create_qphixj_raw_##P##_##T(); \
  QPHIXJ_##P##3_extract_##T##_to_raw((MILCFLOAT *)raw, qphixj, milc2qphixj_parity(parity)); \
  unload_qphixj_raw_##P##_##T##_to_site(dest, raw, parity); \
  destroy_qphixj_raw_##P##_##T(raw); raw = NULL; \
}

/* Map color vector from QPHIXJ field to MILC field */

#define make_unload_to_field(P, T, TYPE, MILCTYPE, MILC_DSTTYPE, MILCFLOAT) \
void \
unload_qphixj_##P##_##T##_to_field( MILC_DSTTYPE *dest, TYPE *qphixj, int parity){ \
  MILCTYPE *raw; \
  raw = create_qphixj_raw_##P##_##T(); \
  QPHIXJ_##P##3_extract_##T##_to_raw((MILCFLOAT *)raw, qphixj, milc2qphixj_parity(parity)); \
  unload_qphixj_raw_##P##_##T##_to_field(dest, raw, parity); \
  destroy_qphixj_raw_##P##_##T(raw); raw = NULL; \
}

#if 0
/* Generate QPHIXJ asqtad links from gauge field in site structure */

#define make_create_L_from_site_gauge(P, TYPE, MILCTYPE, MILCFLOAT) \
TYPE * \
create_qphixj_##P##_L_from_site_gauge( QPHIXJ_info_t *info, \
    QPHIXJ_asqtad_coeffs_t *coeffs, field_offset src, int parity) \
{ \
  MILCTYPE *raw; \
  TYPE *qphixj; \
  QPHIXJ_##P##3_GaugeField *gauge; \
  raw = create_qphixj_raw4_##P##_G_from_site(src, parity); \
  if(raw == NULL)terminate(1); \
  gauge = QPHIXJ_##P##3_create_G_from_raw((MILCFLOAT *)raw, \
				  milc2qphixj_parity(parity)); \
  destroy_qphixj_raw4_##P##_G(raw); raw = NULL; \
  qphixj = QPHIXJ_##P##3_asqtad_create_L_from_G(info, coeffs, gauge); \
  QPHIXJ_##P##3_destroy_G(gauge); gauge = NULL; \
  return qphixj; \
}
#endif

/* Map site gauge links and raw clover term to QPHIXJ field */

#define make_create_L_from_sites(P, TYPE, MILCTYPE, MILCFLOAT) \
TYPE * \
 create_qphixj_##P##_L_from_sites( MILCFLOAT *raw_clov,int parity) \
{ \
  MILCTYPE *raw_links; \
  TYPE *qphixj; \
  raw_links = create_qphixj_raw4_##P##_G_from_site(F_OFFSET(link[0]),parity); \
  qphixj = QPHIXJ_##P##3_wilson_create_L_from_raw((MILCFLOAT *)raw_links, \
	  raw_clov, milc2qphixj_parity(parity)); \
  destroy_G(raw_links); \
  return qphixj; \
}

/* Map site gauge links, raw backward links, and raw clover term to QPHIXJ field */

#define make_create_L_from_fieldback_and_sites(P, TYPE, MILCTYPE, MILCFLOAT) \
TYPE * \
create_qphixj_##P##_L_from_fieldback_and_sites( MILCFLOAT *raw_clov, \
						MILCTYPE *fieldback, int parity) \
{ \
  MILCTYPE *raw_links; \
  TYPE *qphixj; \
  raw_links = create_qphixj_raw4_##P##_G_from_site(F_OFFSET(link[0]),EVENANDODD); \
  qphixj = QPHIXJ_##P##3_wilson_create_L_from_raw((MILCFLOAT *)raw_links, \
	  (MILCFLOAT *)fieldback, raw_clov, milc2qphixj_parity(parity)); \
  destroy_qphixj_raw4_##P##_G(raw_links); \
  return qphixj; \
}

/* Map preconstructed fat and long gauge fields from MILC field to QPHIXJ field */
/* THIS IS UNFINISHED !!! DO NOT USE */

#define make_create_L_from_fields(P, TYPE, MILCTYPE, MILC_SRC_TYPE, MILCFLOAT) \
TYPE * \
 create_qphixj_##P##_L_from_fields( MILC_SRC_TYPE *fat, MILC_SRC_TYPE *lng, \\
    MILC_SRC_TYPE *fatback, MILC_SRC_TYPE *lngback, int parity) \
{ \
  MILCTYPE *rawfat; \
  MILCTYPE *rawlng; \
  MILCTYPE *rawfatback; \
  MILCTYPE *rawlngback; \
  TYPE *qphixj; \
  rawfat = create_qphixj_raw4_##P##_G_from_field(fat, parity); \
  if(rawfat == NULL)terminate(1); \
  rawlng = create_qphixj_raw4_##P##_G_from_field(lng, parity); \
  if(rawlng == NULL)terminate(1); \
  rawfatback = create_qphixj_raw4_##P##_G_from_field(fatback, parity); \
  if(rawfatback == NULL)terminate(1); \
  rawlngback = create_qphixj_raw4_##P##_G_from_field(lngback, parity); \
  if(rawlngback == NULL)terminate(1); \
  qphixj = QPHIXJ_##P##3_asqtad_create_L_from_raw((MILCFLOAT *)rawfat, \
    (MILCFLOAT *)rawlng, (MILCFLOAT *)rawfatback, (MILCFLOAT *)rawlngback, \
    milc2qphixj_parity(parity)); \
  destroy_qphixj_raw4_##P##_G(rawfat); rawfat = NULL; \
  destroy_qphixj_raw4_##P##_G(rawlng); rawlng = NULL; \
  destroy_qphixj_raw4_##P##_G(rawfatback); rawfatback = NULL; \
  destroy_qphixj_raw4_##P##_G(rawlngback); rawlngback = NULL; \
  return qphixj; \
}

/* Extract MILC fat link and long link fields from QPHIXJ asqtad links */
/* THIS IS UNFINISHED !!! DO NOT USE */

#define make_unload_L_to_fields(P, TYPE, MILCTYPE, MILC_DST_TYPE, MILCFLOAT) \
void \
unload_qphixj_##P##_L_to_fields( MILC_DST_TYPE *fat, MILC_DST_TYPE *lng, TYPE* qphixj, \
   int parity) \
{ \
  MILCTYPE *rawfat; \
  MILCTYPE *rawlng; \
  rawfat  = create_qphixj_raw4_##P##_G(); \
  if(rawfat == NULL)terminate(1); \
  rawlng = create_qphixj_raw4_##P##_G(); \
  if(rawlng == NULL)terminate(1); \
  QPHIXJ_##P##3_asqtad_extract_L_to_raw((MILCFLOAT *)rawfat, \
          (MILCFLOAT *)rawlng, qphixj, milc2qphixj_parity(parity)); \
  unload_qphixj_raw4_##P##_G_to_field(fat, rawfat, parity); \
  if(lng != NULL) \
  unload_qphixj_raw4_##P##_G_to_field(lng, rawlng, parity); \
  destroy_qphixj_raw4_##P##_G(rawfat); rawfat = NULL; \
  destroy_qphixj_raw4_##P##_G(rawlng); rawlng = NULL; \
  return; \
}


/* Extract MILC fat link and long link fields from QPHIXJ HISQ links */

#define make_unload_hisq_L_to_fields(P, TYPE, MILCTYPE, MILC_DST_TYPE, MILCFLOAT) \
void \
unload_qphixj_##P##_hisq_L_to_fields( MILC_DST_TYPE *fat, MILC_DST_TYPE *lng, TYPE* qphixj, \
   int parity) \
{ \
  MILCTYPE *rawfat; \
  MILCTYPE *rawlng; \
  rawfat  = create_qphixj_raw4_##P##_G(); \
  if(rawfat == NULL)terminate(1); \
  rawlng = create_qphixj_raw4_##P##_G(); \
  if(rawlng == NULL)terminate(1); \
  QPHIXJ_##P##3_hisq_extract_L_to_raw((MILCFLOAT *)rawfat, \
          (MILCFLOAT *)rawlng, qphixj, milc2qphixj_parity(parity)); \
  unload_qphixj_raw4_##P##_G_to_field(fat, rawfat, parity); \
  if(lng != NULL) \
  unload_qphixj_raw4_##P##_G_to_field(lng, rawlng, parity); \
  destroy_qphixj_raw4_##P##_G(rawfat); rawfat = NULL; \
  destroy_qphixj_raw4_##P##_G(rawlng); rawlng = NULL; \
  return; \
}


/* Map MILC clover term to raw QPHIXJ clover */
/* UNTESTED! */

#define DI(x) (-(1 - (x))/2.)
#define TR(x) ((x)/2.)

#define make_map_milc_clov_to_raw(P, MILCFLOAT) \
void map_milc_clov_to_qphixj_raw_##P(MILCFLOAT *raw_clov, clover *milc_clov ){\
  int i,j,c,status=0;  \
  site *s;\
  MILCFLOAT *r;\
\
  check_layout_compatibility(); \
  create_hash_##T(milc_parity);					\
  FORALLSITES_OMP(i,s,private(j,c,r)){	\
\
    j = hash_m2q_##T[i];	   \
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
  } END_LOOP_OMP \
  if(status==1)terminate(1); \
}

/* map_milc_to_qphixj_all.c */

