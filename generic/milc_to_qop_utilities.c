/******************* milc_to_qop_utilities.c ************************/
/* Functions for mapping MILC data layouts to raw QOP layouts       */
/* C. DeTar 10/19/2005                                              */

#include "generic_includes.h"
#include <qop.h>
#include "../include/generic_qop.h"
#include <string.h>

static int is_qop_initialized = 0;


/* The MILC layout routines list the coordinates and assume 4D */

static int milc_node_number(const int coords[]){
  return node_number(coords[0],coords[1],coords[2],coords[3]);
}

static int milc_node_index(const int coords[]){
  return node_index(coords[0],coords[1],coords[2],coords[3]);
}

QOP_evenodd_t milc2qop_parity(int milc_parity){
  switch(milc_parity){
  case(EVEN):       return QOP_EVEN;
  case(ODD ):       return QOP_ODD;
  case(EVENANDODD): return QOP_EVENODD;
  default:
    printf("milc2qop_parity: Bad MILC parity %d\n", milc_parity);
    terminate(1);
  }
  return (QOP_evenodd_t)-999;
}

int qop2milc_parity(QOP_evenodd_t qop_parity){
  switch(qop_parity){
  case(QOP_EVEN):     return EVEN;
  case(QOP_ODD ):     return ODD;
  case(QOP_EVENODD ): return EVENANDODD;
  default:
    printf("qop2milc_parity: Bad QOP parity %d\n", qop_parity);
    terminate(1);
  }
  return -999;
}

/* Initialize QOP */

QOP_status_t initialize_qop(){
  QOP_status_t status;
  static int latsize[4];
  static QOP_layout_t layout;
  
  latsize[0] = nx;
  latsize[1] = ny;
  latsize[2] = nz;
  latsize[3] = nt;

  if(is_qop_initialized)return QOP_SUCCESS;

  layout.node_number = milc_node_number;
  layout.node_index = milc_node_index;
  layout.latdim = 4;
  layout.latsize = latsize;
  layout.machdim = 4;
  layout.machsize = (int *)get_logical_dimensions();
  layout.this_node = this_node;
  layout.sites_on_node = sites_on_node;

  status = QOP_init(&layout);

  if(status == QOP_SUCCESS)
      is_qop_initialized = 1;

  return status;
}

/* milc_to_qop_utilities */
