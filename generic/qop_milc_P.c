/********************  qop_milc_P.c ***************************************/
/* Implementation of generic QOP API for MILC testing                   */
/* C. DeTar 10/19/2005                                                  */

#if ( QOP_Precision == 1 )
#define MYREAL float
#define MYSU3_VECTOR fsu3_vector
#define MYSU3_MATRIX fsu3_matrix
#define MYWILSON_VECTOR fwilson_vector
#else
#define MYREAL double
#define MYSU3_VECTOR dsu3_vector
#define MYSU3_MATRIX dsu3_matrix
#define MYWILSON_VECTOR dwilson_vector
#endif

#include "generic_includes.h"
#include "../include/generic_qop.h"
#include "../include/generic_qopmilc.h"
#include "../include/qop_milc.h"
#include <string.h>
#include <stdlib.h>

extern int qop_is_initialized;

QOP_ColorVector *QOP_create_V_from_raw(MYREAL *src, QOP_evenodd_t qop_parity){
  MYSU3_VECTOR *v;
  char myname[] = "QOP_create_V_from_raw";
  QOP_ColorVector *qopv;

  if(!qop_is_initialized){
    printf("%s: QOP is not initialized\n",myname);
    return NULL;
  }

  v = (MYSU3_VECTOR *)malloc(sites_on_node*sizeof(MYSU3_VECTOR));
  if(v == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }

  memcpy(v, src, sites_on_node*sizeof(MYSU3_VECTOR));
  
  qopv = (QOP_ColorVector *)malloc(sizeof(QOP_ColorVector));
  if(qopv == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }
  
  qopv->v = v;
  qopv->evenodd = qop_parity;
  return qopv;
}

QOP_DiracFermion *QOP_create_D_from_raw(MYREAL *src, QOP_evenodd_t qop_parity){
  MYWILSON_VECTOR *d;
  char myname[] = "QOP_create_D_from_raw";
  QOP_DiracFermion *qopd;

  if(!qop_is_initialized){
    printf("%s: QOP is not initialized\n",myname);
    return NULL;
  }

  d = (MYWILSON_VECTOR *)malloc(sites_on_node*sizeof(MYWILSON_VECTOR));
  if(d == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }
  
  memcpy(d, src, sites_on_node*sizeof(MYWILSON_VECTOR));

  qopd = (QOP_DiracFermion *)malloc(sizeof(QOP_DiracFermion));
  if(qopd == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }
  
  qopd->d = d;
  qopd->evenodd = qop_parity;
  return qopd;
}


QOP_GaugeField *QOP_create_G_from_raw(MYREAL *src[], QOP_evenodd_t qop_parity){
  MYSU3_MATRIX *g;
  char myname[] = "QOP_create_G_from_raw";
  QOP_GaugeField *qopg;
  int dir;

  if(!qop_is_initialized){
    printf("%s: QOP is not initialized\n",myname);
    return NULL;
  }

  g = (MYSU3_MATRIX *)malloc(sites_on_node*sizeof(MYSU3_MATRIX)*4);
  if(g == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }
  
  FORALLUPDIR(dir){
    memcpy(g+dir*sites_on_node, src[dir], sites_on_node*sizeof(MYSU3_MATRIX));
  }

  qopg = (QOP_GaugeField *)malloc(sizeof(QOP_GaugeField));
  if(qopg == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }
  
  qopg->g = g;
  qopg->evenodd = qop_parity;
  return qopg;
}


QOP_Force *QOP_create_F_from_raw(MYREAL *src[], QOP_evenodd_t qop_parity){
  MYSU3_MATRIX *f;
  char myname[] = "QOP_create_F_from_raw";
  QOP_Force *qopf;
  int dir;

  if(!qop_is_initialized){
    printf("%s: QOP is not initialized\n",myname);
    return NULL;
  }

  f = (MYSU3_MATRIX *)malloc(sites_on_node*sizeof(MYSU3_MATRIX)*4);
  if(f == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }
  
  FORALLUPDIR(dir){
    memcpy(f+dir*sites_on_node, src[dir], sites_on_node*sizeof(MYSU3_MATRIX));
  }

  qopf = (QOP_Force *)malloc(sizeof(QOP_Force));
  if(qopf == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }
  
  qopf->f = f;
  qopf->evenodd = qop_parity;
  return qopf;
}

QOP_FermionLinksAsqtad *QOP_asqtad_create_L_from_raw(MYREAL *fatlinks[], 
				     MYREAL *longlinks[], QOP_evenodd_t qop_parity){
  char myname[] = "QOP_create_L_from_raw";
  QOP_FermionLinksAsqtad *qopl;

  if(!qop_is_initialized){
    printf("%s: QOP is not initialized\n",myname);
    return NULL;
  }

  qopl = (QOP_FermionLinksAsqtad *)malloc(sizeof(QOP_FermionLinksAsqtad));
  if(qopl == NULL){
    printf("%s: No room\n",myname);
    return NULL;
  }
  
  qopl->fat  = QOP_create_G_from_raw((MYREAL **)fatlinks, qop_parity);
  if(qopl->fat == NULL)return NULL;
  qopl->lng = QOP_create_G_from_raw((MYREAL **)longlinks, qop_parity);
  if(qopl->lng == NULL)return NULL;

  qopl->evenodd = qop_parity;
  return qopl;
}


void QOP_extract_V_to_raw(MYREAL *dest, QOP_ColorVector *src, QOP_evenodd_t qop_parity){
  memcpy(dest, src->v, sites_on_node*sizeof(MYSU3_VECTOR));
}

void QOP_extract_D_to_raw(MYREAL *dest, QOP_DiracFermion *src, QOP_evenodd_t qop_parity){
  memcpy(dest, src->d, sites_on_node*sizeof(MYWILSON_VECTOR));
}

void QOP_extract_G_to_raw(MYREAL *dest[], QOP_GaugeField *src, QOP_evenodd_t qop_parity){
  int dir;
  FORALLUPDIR(dir){
    memcpy(dest[dir], src->g+dir*sites_on_node, 
	   sites_on_node*sizeof(MYSU3_MATRIX));
  }
}

void QOP_extract_F_to_raw(MYREAL *dest[], QOP_Force *src, QOP_evenodd_t qop_parity){
  int dir;
  FORALLUPDIR(dir){
    memcpy(dest[dir], src->f+dir*sites_on_node, 
	   sites_on_node*sizeof(MYSU3_MATRIX));
  }
}

void QOP_asqtad_extract_L_to_raw(MYREAL *fatlinks[], MYREAL *longlinks[],
				 QOP_FermionLinksAsqtad *src,
				 QOP_evenodd_t evenodd){
  QOP_extract_G_to_raw(fatlinks, src->fat, evenodd);
  QOP_extract_G_to_raw(longlinks, src->lng, evenodd);
}

void QOP_asqtad_convert_L_to_raw(MYREAL ***fatlinks, MYREAL ***longlinks,
				 QOP_FermionLinksAsqtad *src,
				 QOP_evenodd_t evenodd){
  *fatlinks = (MYREAL **)create_raw4_G();
  *longlinks = (MYREAL **)create_raw4_G();
  QOP_extract_G_to_raw(*fatlinks, src->fat, evenodd);
  QOP_extract_G_to_raw(*longlinks, src->lng, evenodd);
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

void QOP_asqtad_destroy_L(QOP_FermionLinksAsqtad *field){
  if(field == NULL)return;
  QOP_destroy_G(field->fat);
  QOP_destroy_G(field->lng);
  free(field);
}
