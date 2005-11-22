/********************  qop_milc.c ***************************************/
/* Implementation of generic QOP API for MILC testing                   */
/* C. DeTar 10/19/2005                                                  */

#include "../include/qop_milc.h"
#include "generic_includes.h"
#include <string.h>
#include <stdlib.h>

static QOP_layout_t user_layout;
static int qop_is_initialized = 0;

/* We just use the user's layout */

int QOP_node_number_raw(int coords[]){
  return user_layout.node_number(coords);
}

int QOP_node_index_raw_M(int coords[]){
  return user_layout.node_index(coords);
}

int QOP_node_index_raw_V(int coords[]){
  return user_layout.node_index(coords);
}

int QOP_node_index_raw_D(int coords[]){
  return user_layout.node_index(coords);
}

int QOP_node_index_raw_G(int coords[]){
  return user_layout.node_index(coords);
}

int QOP_node_index_raw_F(int coords[]){
  return user_layout.node_index(coords);
}

QOP_status_t QOP_init(QOP_layout_t *layout){

  /* Just copy in the user's layout */
  user_layout = *layout;

  if(user_layout.latdim != 4){
    printf("QOP_init: Can't do %d lattice dimensions\n",
	   user_layout.latdim);
    return QOP_FAIL;
  }
  qop_is_initialized = 1;
  return QOP_SUCCESS;
}

QOP_status_t QOP_finalize(void){
  qop_is_initialized = 0;
  return QOP_SUCCESS;
}

QOP_ColorVector *QOP_create_V_from_raw(Real *src){
  su3_vector *v;
  char myname[] = __func__;
  QOP_ColorVector *qopv;

  if(!qop_is_initialized){
    printf("%s: QOP is not initialized\n",myname);
    return NULL;
  }

  v = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  if(v == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }

  memcpy(v, src, sites_on_node*sizeof(su3_vector));
  
  qopv = (QOP_ColorVector *)malloc(sizeof(QOP_ColorVector));
  if(qopv == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }
  
  qopv->v = v;
  return qopv;
}

QOP_DiracFermion *QOP_create_D_from_raw(Real *src){
  wilson_vector *d;
  char myname[] = __func__;
  QOP_DiracFermion *qopd;

  if(!qop_is_initialized){
    printf("%s: QOP is not initialized\n",myname);
    return NULL;
  }

  d = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector));
  if(d == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }
  
  memcpy(d, src, sites_on_node*sizeof(wilson_vector));

  qopd = (QOP_DiracFermion *)malloc(sizeof(QOP_DiracFermion));
  if(qopd == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }
  
  qopd->d = d;
  return qopd;
}


QOP_GaugeField *QOP_create_G_from_raw(Real *src[]){
  su3_matrix *g;
  char myname[] = __func__;
  QOP_GaugeField *qopg;
  int dir;

  if(!qop_is_initialized){
    printf("%s: QOP is not initialized\n",myname);
    return NULL;
  }

  g = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix)*4);
  if(g == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }
  
  FORALLUPDIR(dir){
    memcpy(g+dir*sites_on_node, src[dir], sites_on_node*sizeof(su3_matrix));
  }

  qopg = (QOP_GaugeField *)malloc(sizeof(QOP_GaugeField));
  if(qopg == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }
  
  qopg->g = g;
  return qopg;
}


QOP_Force *QOP_create_F_from_raw(Real *src[]){
  su3_matrix *f;
  char myname[] = __func__;
  QOP_Force *qopf;
  int dir;

  if(!qop_is_initialized){
    printf("%s: QOP is not initialized\n",myname);
    return NULL;
  }

  f = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix)*4);
  if(f == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }
  
  FORALLUPDIR(dir){
    memcpy(f+dir*sites_on_node, src[dir], sites_on_node*sizeof(su3_matrix));
  }

  qopf = (QOP_Force *)malloc(sizeof(QOP_Force));
  if(qopf == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }
  
  qopf->f = f;
  return qopf;
}

void QOP_extract_V_to_raw(Real *dest, QOP_ColorVector *src){
  memcpy(dest, src->v, sites_on_node*sizeof(su3_vector));
}

void QOP_extract_D_to_raw(Real *dest, QOP_DiracFermion *src){
  memcpy(dest, src->d, sites_on_node*sizeof(wilson_vector));
}

void QOP_extract_G_to_raw(Real *dest[], QOP_GaugeField *src){
  int dir;
  FORALLUPDIR(dir){
    memcpy(dest[dir], src->g+dir*sites_on_node, 
	   sites_on_node*sizeof(su3_matrix));
  }
}

void QOP_extract_F_to_raw(Real *dest[], QOP_Force *src){
  int dir;
  FORALLUPDIR(dir){
    memcpy(dest[dir], src->f+dir*sites_on_node, 
	   sites_on_node*sizeof(su3_matrix));
  }
}

void QOP_destroy_V(QOP_ColorVector *field){
  if(field == NULL)return;
  if(field->v != NULL)free(field->v);
  free(field);
}

void QOP_destroy_D(QOP_DiracFermion *field){
  if(field == NULL)return;
  if(field->d != NULL)free(field->d);
  free(field);
}

void QOP_destroy_G(QOP_GaugeField *field){
  if(field == NULL)return;
  if(field->g != NULL)free(field->g);
  free(field);
}

void QOP_destroy_F(QOP_Force *field){
  if(field == NULL)return;
  if(field->f != NULL)free(field->f);
  free(field);
}

