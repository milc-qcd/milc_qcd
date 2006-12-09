/****************** fermion_links_asqtad_qop.c ***********************/
/* MIMD version 7 */

/* This is the MILC wrapper for SciDAC Level 3 QOP link smearing */
/* These are generic entry points taking the prevailing MILC precision */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/generic_qop.h"
#include "../include/generic_ks_qop.h"


/*********************************************************************/
/* Create fat and long links and qop_links                           */
/*********************************************************************/
/* Wrappers for MILC call to QOP */
void 
load_fn_links( void ){

  if(PRECISION == 1)
    load_fn_links_F();
  else
    load_fn_links_D();
}

#ifdef DM_DU0
/* Wrappers for MILC call to QOP */
void load_fn_links_dmdu0( void ){

  if(PRECISION == 1)
    load_fn_links_dmdu0_F();
  else
    load_fn_links_dmdu0_D();
}
#endif

void
invalidate_fn_links( void )
{
  invalidate_fn_links_F();
  invalidate_fn_links_D();
}
