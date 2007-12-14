/****** fermion_force_asqtad.c  -- ******************/
/* MIMD version 7 */

/* External entry points 

    eo_fermion_force_multi
    ks_multiff_opt_chr

*/

/* Requires these externals:
   
    fermion_force_asqtad_block
    fermion_force_fn_block
    fermion_force_fn_multi_reverse
    fermion_force_fn_multi
*/   

#include "generic_ks_includes.h"	/* definitions files and prototypes */

/* Set default if undeclared */
#ifndef KS_MULTIFF
#define KS_MULTIFF FNMAT
#endif

/**********************************************************************/
/*   Wrapper for fermion force routines with multiple sources         */
/**********************************************************************/
void eo_fermion_force_multi( Real eps, Real *residues, su3_vector **xxx, 
			     int nterms, int prec, ferm_links_t *fn,
			     ks_action_paths *ap ) {
  switch(KS_MULTIFF){
  case ASVEC:
    fermion_force_asqtad_block( eps, residues, xxx, nterms, VECLENGTH, 
				prec, fn, ap );
    break;
  case FNMATREV:
    fermion_force_fn_multi_reverse( eps, residues, xxx, nterms, fn, ap );
    break;
  case FNMAT:
    fermion_force_fn_multi( eps, residues, xxx, nterms, prec, fn, ap );
    break;
  default:
    fermion_force_fn_multi( eps, residues, xxx, nterms, prec, fn, ap );
  }
}

/**********************************************************************/
/*   Accessor for string describing the option                        */
/**********************************************************************/
const char *ks_multiff_opt_chr( void )
{
  switch(KS_MULTIFF){
  case ASVEC:
    return "ASVEC";
    break;
  case FNMATREV:
    return "FNMATREV";
    break;
  case FNMAT:
    return "FNMAT";
    break;
  default:
    return "FNMAT";
  }
  return NULL;
}


