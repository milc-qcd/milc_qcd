/******************* milc_to_quda_utilities.c ************************/
/* For the QUDA/GPU interface */

#include "generic_includes.h"
#include "../include/generic_quda.h"
#include <string.h>

static int is_quda_initialized = 0;


int initialize_quda(void){

  QudaLayout_t layout;
  const int dim[4] = {nx, ny, nz, nt};
  int status = 0;

  if(is_quda_initialized)return status;

  layout.latsize = dim;
  layout.machsize = get_logical_dimensions();
  qudaInit(layout);

  if(status == 0)
    is_quda_initialized = 1;

  return status;

} /* milc_to_quda_utilities */
