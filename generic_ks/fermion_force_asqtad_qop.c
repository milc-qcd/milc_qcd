/******* fermion_force_asqtad_qop.c ****************/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP fermion force routine */
/* 11/2005 C. DeTar and D. Renner */
/* 5/9/07  C. DeTar support precision selection */

/* External entry points in this file

   eo_fermion_force_oneterm
   eo_fermion_force_twoterms
   eo_fermion_force_multi
   ks_multiff_opt_chr

 */

/* Compile with fermion_force_asqtad_qop_[FD].c */

/*
 * $Log: fermion_force_asqtad_qop.c,v $
 * Revision 1.28  2011/11/29 20:45:56  detar
 * Support new fermion links scheme
 *
 * Revision 1.27  2007/12/14 04:36:31  detar
 * Major modification to support HISQ.
 *
 * Revision 1.26  2007/11/09 16:42:41  detar
 * Pull FN link calculation out of inverters
 *
 * Revision 1.25  2007/10/09 20:10:14  detar
 * Add ferm_links_t and ks_action_paths structures and pass them as params
 *
 * Revision 1.24  2007/05/21 05:06:18  detar
 * Support precision selection for fermion force.  Systematize calling.
 *
 * Revision 1.23  2007/03/27 20:59:25  detar
 * Support qopqdp-0.7.7 and later.  Take timing from qopqdp.
 *
 * Revision 1.22  2006/12/16 13:48:45  detar
 * Mixed precision support in QOP MILC.  Support QOP QDP FNMAT
 *
 * Revision 1.21  2006/12/15 02:58:10  detar
 * Make reporting of remapping time optional through a REMAP macro.
 *
 * Revision 1.20  2006/12/15 02:46:16  detar
 * Support Ludmila's addition of Doug's FNMAT fermion force to QOP
 *
 * Revision 1.19  2006/12/12 18:07:16  detar
 * Correct mixed precision features.  Add 1sum variant of the QDP inverter.
 *
 * Revision 1.18  2006/12/09 13:52:39  detar
 * Add mixed precision capability for KS inverter in QOP and QDP
 *
 * Revision 1.17  2006/11/16 04:17:53  detar
 * Fix printf -> node0_printf to avoid duplicate messages.
 *
 * Revision 1.16  2006/11/13 03:05:26  detar
 * Add timing for remapping and make separate from timing for computation.
 *
 * Revision 1.15  2006/11/04 23:41:17  detar
 * Add QOP and QDP support for FN fermion links
 * Create QDP version of fermion_links_fn_multi
 * Add nrestart parameter for ks_congrad
 *
 * Revision 1.14  2006/10/12 03:45:16  detar
 * Move load_qop_asqtad_coeffs to (new) load_qop_asqtad_coeffs.c to
 * prepare for QOP link fattening.
 *
 * Revision 1.13  2006/10/09 03:51:06  detar
 * Collect multi-source fermion force routines in a single file.
 * Move procedures from ks_imp_rhmc/path_transport_field to path_transport.c
 *
 * Revision 1.12  2006/08/29 15:07:08  detar
 * Correct the weights for multi fermion force.  Prepare for QOPQDP.
 *
 * Revision 1.11  2006/08/26 15:35:35  detar
 * Fix nterms assertion.
 *
 * Revision 1.10  2006/08/13 15:02:32  detar
 * Realign procedures to accommodate ks_imp_rhmc code
 * Add Level 3 wrappers and MILC dummy Level 3 implementation for multiple source
 *
 * Revision 1.9  2006/03/11 04:24:22  detar
 * Change to conform to current Level 3 interface for QOP_asqtad_force_multi
 *
 * Revision 1.8  2006/01/30 20:29:14  detar
 * Fix memory leak
 *
 * Revision 1.7  2006/01/28 17:59:08  detar
 * Add FFTIME reporting
 *
 * Revision 1.6  2005/12/09 17:04:01  detar
 * Fix Header entry
 *
 * Revision 1.5  2005/12/09 17:03:16  detar
 * *** empty log message ***
 *
 * Revision 1.4  2005/12/09 17:02:28  detar
 * Add CVS log and header
 *
 */

#include "generic_ks_includes.h"
#include "../include/generic_qop.h"
#include "../include/generic_ks_qop.h"
#include "../include/fermion_links.h"

/* Set default if undeclared */
#ifndef KS_MULTIFF
#define KS_MULTIFF FNMAT
#endif

//static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_ks/fermion_force_asqtad_qop.c,v 1.28 2011/11/29 20:45:56 detar Exp $";

/**********************************************************************/
/* Standard MILC interface for the single-species Asqtad fermion force
   routine */
/**********************************************************************/
void eo_fermion_force_oneterm( Real eps, Real weight, su3_vector *x_off,
			       int prec, fermion_links_t *fl)
{

  if(prec == 1)
    eo_fermion_force_oneterm_F( eps, weight, x_off, fl);
  else
    eo_fermion_force_oneterm_D( eps, weight, x_off, fl);

}

void eo_fermion_force_oneterm_site( Real eps, Real weight, field_offset x_off_site,
				    int prec, fermion_links_t *fl)
{
  su3_vector *x_off = create_v_field_from_site_member(x_off_site);

  eo_fermion_force_oneterm(eps, weight, x_off, prec, fl);

  destroy_v_field(x_off);

}

/**********************************************************************/
/* Standard MILC interface for the two-species Asqtad fermion force
   routine */
/**********************************************************************/
void eo_fermion_force_twoterms( Real eps, Real weight1, Real weight2, 
				su3_vector *x1_off, su3_vector *x2_off,
				int prec, fermion_links_t *fl) 
{

  if(prec == 1)
    eo_fermion_force_twoterms_F( eps, weight1, weight2, x1_off, x2_off, 
				 fl );
  else
    eo_fermion_force_twoterms_D( eps, weight1, weight2, x1_off, x2_off,
				 fl );
}

void eo_fermion_force_twoterms_site( Real eps, Real weight1, Real weight2, 
				     field_offset x1_off_site, field_offset x2_off_site,
				     int prec, fermion_links_t *fl)
{

  su3_vector *x1_off = create_v_field_from_site_member(x1_off_site);
  su3_vector *x2_off = create_v_field_from_site_member(x2_off_site);

  eo_fermion_force_twoterms(eps, weight1, weight2, x1_off, x2_off, prec, fl);

  destroy_v_field(x2_off);
  destroy_v_field(x1_off);

}


/**********************************************************************/
/*   Version for asqtad.  Parallel transport nterms source vectors    */
/**********************************************************************/

static void 
fermion_force_multi( Real eps, Real *residues, 
		     su3_vector **xxx, int nterms, int prec,
		     fermion_links_t *fl ) 
{

  if(prec == 1)
    fermion_force_multi_F( eps, residues, xxx, nterms, fl );
  else
    fermion_force_multi_D( eps, residues, xxx, nterms, fl );

}

/**********************************************************************/
/*   Parallel transport vectors in blocks of veclength. Asqtad only   */
/**********************************************************************/
/* Requires the xxx1 and xxx2 terms in the site structure */

void 
fermion_force_block( Real eps, Real *residues, 
		     su3_vector **xxx, int nterms, int veclength, 
		     int prec, fermion_links_t *fl ) 
{
  if(prec == 1)
    fermion_force_block_F( eps, residues, xxx, nterms, veclength, fl );
  else
    fermion_force_block_D( eps, residues, xxx, nterms, veclength, fl );

}

/**********************************************************************/
/*   Standard MILC interface for fermion force with multiple sources  */
/**********************************************************************/
void eo_fermion_force_multi( Real eps, Real *residues, su3_vector **xxx, 
			     int nterms, int prec, fermion_links_t *fl ) {

  int veclength;
#ifdef VECLENGTH
  veclength = VECLENGTH;
#else
  veclength = 4;
#endif
  /* Set default method to FNMAT. Set vector length for ASVEC, if
     requested */
  QOP_opt_t qop_ff_opt[2] = {
    {.tag = "fnmat_src_min",.value=4},
    {.tag = "veclength",.value=veclength}
  };

  switch(KS_MULTIFF){
  case ASVEC:
    qop_ff_opt[0].value = nterms + 1;  /* set high threshold for FNMAT */
    if(QOP_asqtad_force_set_opts(qop_ff_opt, 2) != QOP_SUCCESS)
      node0_printf("eo_fermion_force_multi: error setting QOP options\n");
    fermion_force_block( eps, residues, xxx, nterms, veclength, prec, fl );
    break;
  default:  /* FNMAT */
    qop_ff_opt[0].value = 4; /* set sensible threshold for FNMAT */
    if(QOP_asqtad_force_set_opts(qop_ff_opt, 2) != QOP_SUCCESS)
      node0_printf("eo_fermion_force_multi: error setting QOP options\n");
    fermion_force_multi( eps, residues, xxx, nterms, prec, fl );
  }
}

/**********************************************************************/
/*   Accessor for string describing the option                        */
/**********************************************************************/
const char *
ks_multiff_opt_chr( void )
{
  switch(KS_MULTIFF){
  case ASVEC:
    return "ASVEC";
    break;
  default:
    return "FNMAT";
  }
  return NULL;
}

