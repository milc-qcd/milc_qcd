/************************ field_utilities.c *****************************/
/* MIMD version 7 */

/* A collection of utilities for creating and copying fields */
/* (Plus a timing utility) */

#include "generic_includes.h"

/*--------------------------------------------------------------------*/
double start_timing(void){
  double dtime;

#ifdef PRTIME
  dtime = -dclock();
#else
  dtime = 0;
#endif
  return dtime;
}

/*--------------------------------------------------------------------*/
void print_timing(double dtime, char *str){

#ifdef PRTIME
  dtime += dclock();
  node0_printf("Time for %s %e\n",str, dtime);  fflush(stdout);
#endif

}

/*--------------------------------------------------------------------*/
#define make_clear_field(ABBREV, T) \
void clear_##ABBREV##_field(T *x){ \
  memset(x,'\0',sites_on_node*sizeof(T)); \
}

#define make_create_field(ABBREV, T) \
T* create_##ABBREV##_field(void){ \
  T *x; \
  x = (T *)malloc(sites_on_node*sizeof(T)); \
  if(x == NULL){ \
    printf("create_field: no room\n"); \
    terminate(1); \
  } \
  clear_##ABBREV##_field(x); \
  return x; \
}

#define make_copy_field(ABBREV, T) \
void copy_##ABBREV##_field(T *dst, T *src){ \
  memcpy(dst, src, sites_on_node*sizeof(T)); \
}

#define make_destroy_field(ABBREV, T) \
void destroy_##ABBREV##_field(T *x){ \
  if(x != NULL) free(x); \
}

#define make_all_field(ABBREV, T) \
  make_clear_field(ABBREV, T); \
  make_create_field(ABBREV, T); \
  make_copy_field(ABBREV, T); \
  make_destroy_field(ABBREV, T);

/*------------------------------------------------------------------*/
/* Create standard create, copy, clear, destroy for standard types  */

make_all_field(c, complex)
make_all_field(v, su3_vector)
make_all_field(m, su3_matrix)
make_all_field(wv, wilson_vector)
make_all_field(swv, spin_wilson_vector)

/*------------------------------------------------------------------*/
/* copy a gauge field in the site structure - an array of four su3_matrices */
void gauge_field_copy(field_offset src, field_offset dest){
  register int i,dir,src2,dest2;
  register site *s;
  FORALLSITES(i,s){
    src2=src; dest2=dest;
    for(dir=XUP;dir<=TUP; dir++){
      su3mat_copy( (su3_matrix *)F_PT(s,src2),
		   (su3_matrix *)F_PT(s,dest2) );
      src2 += sizeof(su3_matrix);
      dest2 += sizeof(su3_matrix);
    }
  }
}

/*------------------------------------------------------------------*/
/* copy a gauge field from the site structure to a field - an array of
   four su3_matrices.  We imitate the raw QOP format here, but use
   the MILC layout site order . */

su3_matrix **gauge_field_copy_site_to_field(field_offset src){
  int i,dir;
  site *s;
  su3_matrix *tmp;
  su3_matrix **dest;

  dest = (su3_matrix **)malloc(4*sizeof(su3_matrix *));
  FORALLUPDIR(dir){
    dest[dir] = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    if(dest[dir] == NULL){
      printf("gauge_field_copy_site_to_field: No room\n");
      return NULL;
    }
  }

  FORALLSITES(i,s){
    tmp = (su3_matrix *)F_PT(s, src);		\
    FORALLUPDIR(dir){
      su3mat_copy( tmp + dir, dest[dir] + i );
    }
  }

  return dest;
}


/*------------------------------------------------------------------*/
/* copy a gauge field from a field to the site structure.
   The opposite of gauge_field_copy_site_to_field. */

void gauge_field_copy_field_to_site(su3_matrix **src, field_offset dest){
  int i,dir;
  site *s;
  su3_matrix *tmp;

  FORALLSITES(i,s){
    tmp = (su3_matrix *)F_PT(s, dest);		\
    FORALLUPDIR(dir){
      su3mat_copy( src[dir] + i, tmp + dir );
    }
  }
}


/*------------------------------------------------------------------*/
/* destroy a gauge field with the format used above                 */

void destroy_gauge_field(su3_matrix **f){
  int dir;

  if(f != NULL){
    FORALLUPDIR(dir){
      if(f[dir] != NULL)
	free(f[dir]);
    }
    free(f);
  }
}

/*--------------------------------------------------------------------*/
su3_vector *create_v_field_from_site_member(field_offset sv){
  su3_vector *v;
  int i; site *s;

  v = create_v_field();

  FORALLSITES(i,s){
    memcpy(v+i, F_PT(s,sv), sizeof(su3_vector));
  }

  return v;
}

/*--------------------------------------------------------------------*/
void copy_site_member_from_v_field(field_offset sv, su3_vector *v){
  int i; site *s;

  FORALLSITES(i,s){
    memcpy(F_PT(s,sv), v+i, sizeof(su3_vector));
  }
}

/*--------------------------------------------------------------------*/
void add_v_fields(su3_vector *vsum, su3_vector *v1, su3_vector *v2){
  int i;
  site *s;

  FORALLSITES(i,s){
    add_su3_vector(v1+i, v2+i, vsum+i);
  }
}

/*--------------------------------------------------------------------*/
void extract_c_from_v(complex *c, su3_vector *v, int color){
  int i; site *s;
  FORALLSITES(i,s){
    c[i] = v[i].c[color];
  }
}

/*--------------------------------------------------------------------*/
/* Insert complex field into su3_vector field */
void insert_v_from_c(su3_vector *v, complex *c, int color){
  int i; site *s;
  FORALLSITES(i,s){
    v[i].c[color] = c[i];
  }
}
/*--------------------------------------------------------------------*/
ks_prop_field *create_ksp_field(int nc){
  ks_prop_field *ksp;
  int color;
  
  ksp = (ks_prop_field *)malloc(sizeof(ks_prop_field));
  if(ksp == NULL){
    printf("create_ksp_field(%d): No room for temporary\n",this_node);
    terminate(1);
  }

  ksp->nc = nc;

  ksp->v = (su3_vector **) malloc(nc*sizeof(su3_vector *));

  if(ksp->v == NULL){
    printf("create_ksp_field(%d): No room for temporary\n",this_node);
    terminate(1);
  }

  for(color= 0; color < nc; color++)
    ksp->v[color] = create_v_field();
  
  return ksp;
}

/*--------------------------------------------------------------------*/
void clear_ksp_field(ks_prop_field *ksp){
  int color;
  
  for(color= 0; color < ksp->nc; color++)
    clear_v_field(ksp->v[color]);
}

/*--------------------------------------------------------------------*/
void copy_ksp_field(ks_prop_field *kspcopy, ks_prop_field *ksp){
  int color;

  kspcopy->nc = ksp->nc;
  for(color = 0; color < ksp->nc; color++){
    copy_v_field(kspcopy->v[color], ksp->v[color]);
  }
}

/*--------------------------------------------------------------------*/
ks_prop_field *create_ksp_field_copy(ks_prop_field *k){
  ks_prop_field *ksp;

  ksp = create_ksp_field(k->nc);
  copy_ksp_field(ksp, k);

  return ksp;
}

/*--------------------------------------------------------------------*/
void destroy_ksp_field(ks_prop_field *ksp){
  int color;

  for(color = 0; color < ksp->nc; color++){
    destroy_v_field(ksp->v[color]);
  }
  free(ksp->v);
}

/*--------------------------------------------------------------------*/
void copy_v_from_ksp(su3_vector *v, ks_prop_field *ksp, int color){
  int i;
  site *s;

  FORALLSITES(i,s){
    v[i] = ksp->v[color][i];
  }
}

/*--------------------------------------------------------------------*/
void insert_ksp_from_v(ks_prop_field *ksp, su3_vector *v, int color){
  int i;
  site *s;

  FORALLSITES(i,s){
    ksp->v[color][i] = v[i];
  }
}

/*--------------------------------------------------------------------*/
/* Copy complex field from wilson_vector field */
void extract_c_from_wv(complex *c, wilson_vector *wv, 
		       int spin, int color){
  int i; site *s;
  FORALLSITES(i,s){
    c[i] = wv[i].d[spin].c[color];
  }
}

/*--------------------------------------------------------------------*/
/* Insert c field into wv field  */
void insert_wv_from_c(wilson_vector *wv, complex *c, 
		      int spin, int color){
  int i; site *s;
  FORALLSITES(i,s){
    wv[i].d[spin].c[color] = c[i];
  }
}

/*--------------------------------------------------------------------*/
/* Insert v field into wv field  */
void insert_wv_from_v(wilson_vector *wv, su3_vector *v, int spin){
  int i; site *s;
  int c;
  FORALLSITES(i,s){
    for(c = 0; c < 3; c++)
      wv[i].d[spin].c[c] = v[i].c[c];
  }
}

/*--------------------------------------------------------------------*/
/* Insert v field into wv field  */
void extract_v_from_wv(su3_vector *v, wilson_vector *wv, int spin){
  int i; site *s;
  int c;
  FORALLSITES(i,s){
    for(c = 0; c < 3; c++)
      v[i].c[c] = wv[i].d[spin].c[c];
  }
}

/*--------------------------------------------------------------------*/
spin_wilson_vector *extract_swv_from_wp(wilson_prop_field *wp, int color){
  return wp->swv[color];
}

/*--------------------------------------------------------------------*/
/* Reallocates space for the spin wilson vectors */
void rebuild_wp_field(wilson_prop_field *wp){

  int color;

  if(wp->nc <= 0){
    node0_printf("rebuild_wp_field: Illegal %d colors\n", wp->nc);
    terminate(1);
  }
  
  if(wp->swv == NULL){
    printf("rebuild_wp_field(%d): Needs allocated swv\n",this_node);
    terminate(1);
  }

  for(color= 0; color < wp->nc; color++)
    wp->swv[color] = create_swv_field();
}

/*--------------------------------------------------------------------*/
wilson_prop_field *create_wp_field(int nc){
  wilson_prop_field *wp= NULL;

  if(nc <= 0){
    node0_printf("create_wp_field: Illegal %d colors\n", nc);
    terminate(1);
  }
  
  wp = (wilson_prop_field *) malloc(sizeof(wilson_prop_field));
  if(wp == NULL){
    printf("create_wp_field(%d): No room\n",this_node);
    terminate(1);
  }

  wp->nc = nc;
  wp->swv = (spin_wilson_vector **) malloc(nc*sizeof(spin_wilson_vector *));
  if(wp->swv == NULL){
    printf("create_wp_field(%d): No room for %d colors in field\n",this_node, nc);
    terminate(1);
  }

  rebuild_wp_field(wp);
  
  return wp;
}

/*--------------------------------------------------------------------*/
void clear_wp_field(wilson_prop_field *wp){
  int color;
  
  for(color= 0; color < wp->nc; color++)
    clear_swv_field(wp->swv[color]);
}

/*--------------------------------------------------------------------*/
void copy_wp_field(wilson_prop_field *wpcopy, wilson_prop_field *wp){
  int color, i;
  site *s;

  wpcopy->nc = wp->nc;
  for(color = 0; color < wp->nc; color++){
    FORALLSITES(i,s){
      wpcopy->swv[color][i] = wp->swv[color][i];
    }
  }
}

/*--------------------------------------------------------------------*/
wilson_prop_field *create_wp_field_copy(wilson_prop_field *w){
  wilson_prop_field *wp;

  wp = create_wp_field(w->nc);
  copy_wp_field(wp, w);

  return wp;
}

/*--------------------------------------------------------------------*/
/* Transpose source and sink color and spin indices in place */
void transpose_wp_field(wilson_prop_field *wp){
  int c1, c2, s1, s2, i;
  site *s;
  spin_wilson_vector swv[3];

  if(wp->nc != 3){
    node0_printf("transpose_wp_field: Requires three colors to do a transpose, but have %d\n",wp->nc);
    terminate(1);
  }
  FORALLSITES(i,s){
    for(c2 = 0; c2 < 3; c2++)
      swv[c2] = wp->swv[c2][i];
    for(c2 = 0; c2 < 3; c2++)
      for(s2 = 0; s2 < 4; s2++)
	for(s1 = 0; s1 < 4; s1++)
	  for(c1 = 0; c1 < 3; c1++)
	    wp->swv[c1][i].d[s1].d[s2].c[c2] = swv[c2].d[s2].d[s1].c[c1];
  }
}

/*--------------------------------------------------------------------*/
void copy_wv_from_v(wilson_vector *wv, su3_vector *v, int spin){
  int i;
  
  FORALLFIELDSITES(i){
    wv[i].d[spin] = v[i];
  }
}

/*--------------------------------------------------------------------*/
void copy_wv_from_wp(wilson_vector *wv, wilson_prop_field *wp, 
		     int color, int spin){
  int i;
  site *s;
  
  FORALLSITES(i,s){
    wv[i] = wp->swv[color][i].d[spin];
  }
}

/*--------------------------------------------------------------------*/
void copy_wv_from_wprop(wilson_vector *wv, wilson_propagator *wprop, 
			int color, int spin){
  int i;
  site *s;
  
  FORALLSITES(i,s){
    wv[i] = wprop[i].c[color].d[spin];
  }
}

/*--------------------------------------------------------------------*/
void copy_wv_from_swv(wilson_vector *wv, spin_wilson_vector *swv, int spin){
  int i;
  site *s;
  
  FORALLSITES(i,s){
    wv[i] = swv[i].d[spin];
  }
}

/*--------------------------------------------------------------------*/
void copy_wp_from_wv(wilson_prop_field *wp, wilson_vector *wv, 
		     int color, int spin){
  int i;
  site *s;
  
  FORALLSITES(i,s){
    wp->swv[color][i].d[spin] = wv[i];
  }
}

/*--------------------------------------------------------------------*/
/* Free major storage, but preserves nc and swv */
void free_wp_field(wilson_prop_field *wp){
  int color;

  if(wp == NULL)return;
  for(color = 0; color < wp->nc; color++){
    destroy_swv_field(wp->swv[color]);
    wp->swv[color] = NULL;
  }
}
/*--------------------------------------------------------------------*/
/* Free everything */
void destroy_wp_field(wilson_prop_field *wp){

  if(wp == NULL)return;
  free_wp_field(wp);
  if(wp->swv != NULL)
    free(wp->swv);
  free(wp);
}

/*--------------------------------------------------------------------*/

wilson_vector *create_wv_from_swv(spin_wilson_vector *swv, int spin){
  
  int i;
  site *s;
  wilson_vector *wv = create_wv_field();
  
  FORALLSITES(i,s){
    wv[i] = swv[i].d[spin];
  }
  return wv;
}

/*--------------------------------------------------------------------*/

void insert_swv_from_wv(spin_wilson_vector *swv, int spin, wilson_vector *wv){
  
  int i;
  site *s;
  
  FORALLSITES(i,s){
    swv[i].d[spin] = wv[i];
  }

}

