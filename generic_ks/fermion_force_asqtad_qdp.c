/****** fermion_force_asqtad_qdp.c  -- ******************/
/* MIMD version 7 */
/* Fermion force for Asqtad optimized with single direction QDP operations

 * C.D. 5/07 Calculate with specified precision
 *           Main code moved to fermion_force_asqtad_qdp_P.c
 */

/* External entry points in this file

   eo_fermion_force_oneterm
   eo_fermion_force_twoterms
   eo_fermion_force_multi
   fermion_force_asqtad_block
   fermion_force_asqtad_multi
   ks_multiff_opt_chr

*/


/* Requires these externals 

   fermion_force_fn_multi

*/

/* Compile with fermion_force_fn_multi_qdp.c and 
   fermion_force_asqtad_qdp_[FD].c */

#include "generic_ks_includes.h"
#include "../include/generic_qdp.h"
#include "../include/generic_ks_qdp.h"

/* Set default if undeclared */
#ifndef KS_MULTIFF
#define KS_MULTIFF FNMAT
#endif

/* Standard MILC interface for the single-species Asqtad fermion force routine */
void eo_fermion_force_oneterm( Real eps, Real weight, field_offset x_off,
			       int prec, ferm_links_t *fn, ks_action_paths *ap )
{
  if( prec == 1 )
    eo_fermion_force_oneterm_F( eps, weight, x_off, fn, ap );
  else
    eo_fermion_force_oneterm_D( eps, weight, x_off, fn, ap );
}

/* Standard MILC interface for the two-species Asqtad fermion force routine */

void eo_fermion_force_twoterms( Real eps, Real weight1, Real weight2, 
				field_offset x1_off, field_offset x2_off,
				int prec, ferm_links_t *fn, 
				ks_action_paths *ap ) 
{
  if( prec == 1 )
    eo_fermion_force_twoterms_F( eps, weight1, weight2, x1_off, x2_off,
				 fn, ap );
  else
    eo_fermion_force_twoterms_D( eps, weight1, weight2, x1_off, x2_off,
				 fn, ap );
}

/**********************************************************************/
/*   Parallel transport vectors in blocks of veclength. Asqtad only   */
/**********************************************************************/
/* Requires the xxx1 and xxx2 terms in the site structure */

void 
fermion_force_asqtad_block( Real eps, Real *residues, 
			    su3_vector **xxx, int nterms, int veclength, 
			    int prec, ferm_links_t *fn, ks_action_paths *ap) 
{
  if( prec == 1)
    fermion_force_asqtad_block_F( eps, residues, xxx, nterms, veclength,
				  fn, ap);
  else
    fermion_force_asqtad_block_D( eps, residues, xxx, nterms, veclength, 
				  fn, ap );
}


/**********************************************************************/
/*   Wrapper for fermion force routines with multiple sources         */
/**********************************************************************/
void 
fermion_force_asqtad_multi( Real eps, Real *residues, 
			    su3_vector **xxx, int nterms, int prec,
			    ferm_links_t *fn, ks_action_paths *ap ) 
{
  if( prec == 1)
    fermion_force_asqtad_multi_F( eps, residues, xxx, nterms, fn, ap );
  else
    fermion_force_asqtad_multi_D( eps, residues, xxx, nterms, fn, ap );
}

/**********************************************************************/
/*   Wrapper for fermion force routines with multiple sources         */
/**********************************************************************/
void eo_fermion_force_multi( Real eps, Real *residues, 
			     su3_vector **xxx, int nterms, int prec,
			     ferm_links_t *fn, ks_action_paths *ap ) 
{
  switch(KS_MULTIFF){
  case ASVEC:
    fermion_force_asqtad_block( eps, residues, xxx, nterms, VECLENGTH,
				prec, fn, ap );
    break;
  case FNMAT:
    fermion_force_fn_multi( eps, residues, xxx, nterms, prec, fn, ap );
    break;
  default:
    fermion_force_asqtad_multi( eps, residues, xxx, nterms, prec, fn, ap );
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
  case FNMAT:
    return "FNMAT";
    break;
  default:
    return "FNMAT (default choice)";
  }
  return NULL;
}

