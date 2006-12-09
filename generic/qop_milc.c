/********************  qop_milc.c ***************************************/
/* Implementation of generic QOP API for MILC testing                   */
/* C. DeTar 10/19/2005                                                  */

#include "generic_includes.h"
#include "../include/qop_milc.h"
#include <string.h>
#include <stdlib.h>

static QOP_layout_t user_layout;
int qop_is_initialized = 0;

/* We just use the user's layout */

int QOP_node_number_raw(int coords[]){
  return user_layout.node_number(coords);
}

int QOP_node_index_raw_M(int coords[], QOP_evenodd_t evenodd){
  return user_layout.node_index(coords);
}

int QOP_node_index_raw_V(int coords[], QOP_evenodd_t evenodd){
  return user_layout.node_index(coords);
}

int QOP_node_index_raw_D(int coords[], QOP_evenodd_t evenodd){
  return user_layout.node_index(coords);
}

int QOP_node_index_raw_G(int coords[], QOP_evenodd_t evenodd){
  return user_layout.node_index(coords);
}

int QOP_node_index_raw_F(int coords[], QOP_evenodd_t evenodd){
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

