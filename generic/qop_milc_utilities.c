/******************* qop_milc_utilities.c ************************/
/* MILC Version 7 */

/* Functions for mapping QOP-MILC data of a specific precision to flat
   MILC arrays of the prevailing precision specified by the PRECISION
   macro.  Supports the MILC test version of QOP in mixed precision
   calculation */

/* C. DeTar 12/15/2006                                              */

#include "generic_includes.h"
#include <string.h>

/* Convert (or copy) MILC types between specific and prevailing precision */

#if (PRECISION == 1)

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

/********************************************************************/
/* su3_matrix field conversion                                      */
/********************************************************************/

#if ( PRECISION == 1 )

/* No copy necessary if the prevailing precision matches the input array */

su3_matrix *
create_links_from_qop_milc_F(fsu3_matrix *src)
{
  return src;
}

void
destroy_links_from_qop_milc_F(su3_matrix *g){
}

su3_matrix *
create_links_from_qop_milc_D(dsu3_matrix *src)
{
  su3_matrix *g;
  int i, dir;
  site *s;

  g = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node*4);
  if(g == NULL){
    printf("create_links_from_qop_milc_D(%d): No room\n",this_node);
    return NULL;
  }
  FORALLUPDIR(dir){
    FORALLSITES(i,s){
      d2p_mat(g+dir*sites_on_node+i, src+dir*sites_on_node+i);
    }
  }
  return g;
}

void
destroy_links_from_qop_milc_D(su3_matrix *g){
  if(g != NULL) free(g);
}


#else

su3_matrix *
create_links_from_qop_milc_F(fsu3_matrix *src)
{
  su3_matrix *g;
  int i, dir;
  site *s;

  g = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node*4);
  if(g == NULL){
    printf("create_links_from_qop_milc_F(%d): No room\n",this_node);
    return NULL;
  }
  FORALLUPDIR(dir){
    FORALLSITES(i,s){
      f2p_mat(g+dir*sites_on_node+i, src+dir*sites_on_node+i);
    }
  }
  return g;
}

void
destroy_links_from_qop_milc_F(su3_matrix *g){
  if(g != NULL) free(g);
}

su3_matrix *
create_links_from_qop_milc_D(dsu3_matrix *src)
{
  return src;
}

void
destroy_links_from_qop_milc_D(su3_matrix *g){
}


#endif

/********************************************************************/
/* su3_vector field conversion                                      */
/********************************************************************/

#if ( PRECISION == 1 )

su3_vector *
create_latvec_from_qop_milc_F(fsu3_vector *src)
{
  return src;
}

void
destroy_latvec_from_qop_milc_F(su3_vector *v){
  return;
}

su3_vector *
create_latvec_from_qop_milc_D(dsu3_vector *src)
{
  su3_vector *v;
  int i;
  site *s;

  v = (su3_vector *)malloc(sizeof(su3_vector)*sites_on_node);
  if(v == NULL){
    printf("create_latvec_from_qop_milc_D(%d): No room\n",this_node);
    return NULL;
  }
  FORALLSITES(i,s){
    d2p_vec(v + i, src + i);
  }
  return v;
}

void
destroy_latvec_from_qop_milc_D(su3_vector *v){
  if(v != NULL) free(v);
}


#else

su3_vector *
create_latvec_from_qop_milc_F(fsu3_vector *src)
{
  su3_vector *v;
  int i;
  site *s;

  v = (su3_vector *)malloc(sizeof(su3_vector)*sites_on_node);
  if(v == NULL){
    printf("create_latvec_from_qop_milc_F(%d): No room\n",this_node);
    return NULL;
  }
  FORALLSITES(i,s){
    f2p_vec(v + i, src + i);
  }
  return v;
}

void
destroy_latvec_from_qop_milc_F(su3_vector *v){
  if(v != NULL) free(v);
}

su3_vector *
create_latvec_from_qop_milc_D(dsu3_vector *src)
{
  return src;
}

void
destroy_latvec_from_qop_milc_D(su3_vector *v){
  return;
}

#endif

/********************************************************************/
/* su3_vector field copy from prevailing MILC precision to specific
   QOP precision */
/********************************************************************/

#if ( PRECISION == 1 )
void
copy_latvec_to_qop_milc_F( fsu3_vector *dest, su3_vector *src)
{
  if(dest != src)
    memcpy(dest, src, sites_on_node*sizeof(fsu3_vector));
}

void
copy_latvec_to_qop_milc_D( dsu3_vector *dest, su3_vector *src)
{
  int i;
  site *s;

  FORALLSITES(i,s){
    p2d_vec(dest + i, src + i);
  }
}

#else

void
copy_latvec_to_qop_milc_F(fsu3_vector *dest, su3_vector *src)
{
  int i;
  site *s;

  FORALLSITES(i,s){
    p2f_vec(dest + i, src + i);
  }
}

void
copy_latvec_to_qop_milc_D(dsu3_vector *dest, su3_vector *src)
{
  if(dest != src)
    memcpy(dest, src, sites_on_node*sizeof(dsu3_vector));
}

#endif

/* qop_milc_utilities */
