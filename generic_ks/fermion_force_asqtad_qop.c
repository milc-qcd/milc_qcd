/******* fermion_force_asqtad_milc_qop.c ****************/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP fermion force routine */
/* 11/2005 C. DeTar and D. Renner */

/*
 * $Log: fermion_force_asqtad_qop.c,v $
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

#if (QOP_Precision==1)
#define LOAD_QOP_ASQTAD_COEFFS load_qop_F_asqtad_coeffs
#else
#define LOAD_QOP_ASQTAD_COEFFS load_qop_D_asqtad_coeffs
#endif

#include "generic_ks_includes.h"
#include "../include/generic_qop.h"
#include "../include/generic_ks_qop.h"

/* Set default if undeclared */
#ifndef KS_MULTIFF
#define KS_MULTIFF FNMAT
#endif

static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_ks/fermion_force_asqtad_qop.c,v 1.21 2006/12/15 02:58:10 detar Exp $";

/**********************************************************************/
/* Standard MILC interface for the single-species Asqtad fermion force
   routine */
/**********************************************************************/
void eo_fermion_force_oneterm( Real eps, Real weight, field_offset x_off )
{

  /* For example weight = nflavors/4 */

  su3_matrix **rawlinks;
  su3_matrix **rawmom;
  su3_vector *rawvecx;

  QOP_GaugeField *links;
  QOP_Force *mom;
  QOP_ColorVector *vecx;
  
  QOP_asqtad_coeffs_t coeff;
  QOP_info_t info;

#ifdef FFTIME
  int nflop = 253935;
#endif
  double dtime;
  double remaptime = -dclock();

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    node0_printf("eo_fermion_force_oneterm: Error initializing QOP\n");
    terminate(1);
  }

  /* Load gauge links and momentum */
  load_links_and_mom_site( &links, &mom, &rawlinks, &rawmom );

  /* Copy color vector from site structure to raw and then to QOP format */
  rawvecx = create_raw_V_from_site(x_off,EVENANDODD);
  if(rawvecx == NULL)terminate(1);
  vecx = QOP_create_V_from_raw((Real *)rawvecx,QOP_EVENODD);

  /* Load coefficients */
  LOAD_QOP_ASQTAD_COEFFS(&coeff, weight, get_quark_path_coeff());

  /* Compute fermion force */
  dtime = -dclock();
  QOP_asqtad_force(&info, links, mom, &coeff, eps, vecx);
  dtime += dclock();
  remaptime -= dtime;

  /* Unload momentum and destroy storage for momentum and links */
  unload_links_and_mom_site(  &links, &mom, &rawlinks, &rawmom );

  /* Free QOP source vector */
  QOP_destroy_V(vecx);  vecx = NULL;

  /* Free raw source vector */
  destroy_raw_V(rawvecx); rawvecx = NULL;

  remaptime += dclock();
#ifdef FFTIME
  node0_printf("FFTIME:  time = %e (qop) terms = 1 mflops = %e\n",dtime,
	       (Real)nflop*volume/(1e6*dtime*numnodes()) );
#ifdef REMAP
  node0_printf("FFREMAP:  time = %e\n",remaptime);
#endif
#endif
}

/**********************************************************************/
/* Standard MILC interface for the two-species Asqtad fermion force
   routine */
/**********************************************************************/
void eo_fermion_force_twoterms( Real eps, Real weight1, Real weight2, 
			   field_offset x1_off, field_offset x2_off ) {

  /* For example weight1 = nflavor1/4; weight2 = nflavor2/4 */
  su3_matrix **rawlinks;
  su3_matrix **rawmom;
  su3_vector *rawvecx[2];

  QOP_GaugeField *links;
  QOP_Force *mom;
  QOP_ColorVector *vecx[2];
  
  QOP_asqtad_coeffs_t coeff;
  Real epsv[2];
  QOP_info_t info;
  QOP_opt_t qop_ff_opt = {.tag = "st"};

#ifdef FFTIME
  int nflop = 433968;
#endif
  double dtime;
  double remaptime = -dclock();

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("eo_fermion_force_twoterms: Error initializing QOP\n");
    terminate(1);
  }

  /* Load gauge links and momentum */
  load_links_and_mom_site( &links, &mom, &rawlinks, &rawmom );

  /* Loop over fermion species */
  /* Copy color vectors from site structure to raw and then to QOP format */
  rawvecx[0] = create_raw_V_from_site(x1_off, EVENANDODD);
  if(rawvecx[0] == NULL)terminate(1);
  vecx[0] = QOP_create_V_from_raw((Real *)rawvecx[0],QOP_EVENODD);

  rawvecx[1] = create_raw_V_from_site(x2_off, EVENANDODD);
  if(rawvecx[1] == NULL)terminate(1);
  vecx[1] = QOP_create_V_from_raw((Real *)rawvecx[1],QOP_EVENODD);

  /* Load coefficients */
  epsv[0] = eps*weight1;  epsv[1] = eps*weight2;
  LOAD_QOP_ASQTAD_COEFFS(&coeff, 1., get_quark_path_coeff());

  /* Compute fermion force */
  dtime = -dclock();
  qop_ff_opt.value = 0;  /* ASVEC method is appropriate for 2 sources */
  QOP_asqtad_force_set_opts(&qop_ff_opt, 1);
  QOP_asqtad_force_multi(&info, links, mom, &coeff, epsv, vecx, 2);
  dtime += dclock();
  remaptime -= dtime;

  /* Unload momentum and destroy storage for momentum and links */
  unload_links_and_mom_site(  &links, &mom, &rawlinks, &rawmom );

  /* Free QOP source vectors */
  QOP_destroy_V(vecx[0]);  vecx[0] = NULL;
  QOP_destroy_V(vecx[1]);  vecx[1] = NULL;

  /* Free raw source vectors */
  destroy_raw_V(rawvecx[0]); rawvecx[0] = NULL;
  destroy_raw_V(rawvecx[1]); rawvecx[1] = NULL;

  remaptime += dclock();
#ifdef FFTIME
  node0_printf("FFTIME:  time = %e (qop) terms = 2 mflops = %e\n",dtime,
	       (Real)nflop*volume/(1e6*dtime*numnodes()) );
#ifdef REMAP
  node0_printf("FFREMAP:  time = %e\n",remaptime);
#endif
#endif
}


/**********************************************************************/
/*   Version for asqtad.  Parallel transport nterms source vectors    */
/**********************************************************************/

void eo_fermion_force_asqtad_multi( Real eps, Real *residues, 
    su3_vector **xxx, int nterms ) {
  su3_matrix **rawlinks;
  su3_matrix **rawmom;
  su3_vector *rawvecx = NULL;

  QOP_GaugeField *links;
  QOP_Force *mom;
  QOP_ColorVector **vecx;
  
  QOP_asqtad_coeffs_t coeff;
  int i;
  Real *epsv;
  QOP_info_t info;
  char myname[] = "eo_fermion_force_asqtad_multi";

#ifdef FFTIME
  int nflop = 433968;
#endif
  double dtime;
  double remaptime = -dclock();

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("%s(%d): Error initializing QOP\n",myname,this_node);
    terminate(1);
  }

  /* Load gauge links and momentum */
  load_links_and_mom_site( &links, &mom, &rawlinks, &rawmom );

  /* Make space for QOP vector pointers */
  vecx = (QOP_ColorVector **)malloc(sizeof(QOP_ColorVector *)*nterms);
  if(vecx == NULL){
    printf("%s(%d): no room for vecx pointers\n",myname,this_node);
    terminate(1);
  }

  /* Copy color vectors from field to raw and then to QOP format */
  for(i = 0; i < nterms; i++){
    rawvecx = create_raw_V_from_field(xxx[i], EVENANDODD);
    if(rawvecx == NULL)terminate(1);
    vecx[i] = QOP_create_V_from_raw((Real *)rawvecx,QOP_EVENODD);
    destroy_raw_V(rawvecx); rawvecx = NULL;
  }

  free(rawvecx);

  /* Make space for weights */
  epsv = (Real *)malloc(sizeof(Real)*nterms);
  /* Load coefficients */
  for(i = 0; i < nterms; i++) epsv[i] = eps*residues[i];
  LOAD_QOP_ASQTAD_COEFFS(&coeff, 1., get_quark_path_coeff());

  /* Compute fermion force */
  dtime = -dclock();
  QOP_asqtad_force_multi(&info, links, mom, &coeff, epsv, vecx, nterms);
  dtime += dclock();
  remaptime -= dtime;

  /* Unload momentum and destroy storage for momentum and links */
  unload_links_and_mom_site(  &links, &mom, &rawlinks, &rawmom );

  /* Free QOP source vectors */
  for(i = 0; i < nterms; i++){
    QOP_destroy_V(vecx[i]);  vecx[i] = NULL;
  }

  free(epsv);

  remaptime += dclock();
#ifdef FFTIME
  node0_printf("FFTIME:  time = %e (qop ASVEC) terms = %d mflops = %e\n",dtime,
	       nterms, (Real)nflop*volume/(1e6*dtime*numnodes()) );
#ifdef REMAP
  node0_printf("FFREMAP:  time = %e\n",remaptime);
#endif
#endif

}

/**********************************************************************/
/*   Parallel transport vectors in blocks of veclength. Asqtad only   */
/**********************************************************************/
/* Requires the xxx1 and xxx2 terms in the site structure */

void eo_fermion_force_asqtad_block( Real eps, Real *residues, 
	    su3_vector **xxx, int nterms, int veclength ) {

  int i,j;
  site *s;

  /* First do blocks of size veclength */
  for( j = 0;  j <= nterms-veclength; j += veclength )
    eo_fermion_force_asqtad_multi( eps, &(residues[j]), xxx+j, veclength );
  
#ifndef ONEMASS
  /* Continue with pairs if needed */
  for( ; j <= nterms-2 ; j+=2 ){
    FORALLSITES(i,s){
      s->xxx1 = xxx[j  ][i] ;
      s->xxx2 = xxx[j+1][i] ;
    }
    eo_fermion_force_twoterms( eps, residues[j], residues[j+1],
			       F_OFFSET(xxx1), F_OFFSET(xxx2) );
  }

  /* Finish with a single if needed */
  for( ; j <= nterms-1; j++ ){
    FORALLSITES(i,s){ s->xxx1 = xxx[j][i] ; }
    eo_fermion_force_oneterm( eps, residues[j], F_OFFSET(xxx1) );
  }
#else
  /* Thrown in just to make it compile.  Probably never used. */
  /* Finish with a single if needed */
  for( ; j <= nterms-1; j++ ){
    FORALLSITES(i,s){ s->xxx = xxx[j][i] ; }
    eo_fermion_force_oneterm( eps, residues[j], F_OFFSET(xxx) );
  }
#endif

}

/**********************************************************************/
/*   Standard MILC interface for fermion force with multiple sources  */
/**********************************************************************/
void eo_fermion_force_multi( Real eps, Real *residues, su3_vector **xxx, 
			     int nterms ) {

  int veclength;
#ifdef VECLENGTH
  veclength = VECLENGTH;
#else
  veclength = 4;
#endif
  QOP_opt_t qop_ff_opt = {.tag = "st"};

  switch(KS_MULTIFF){
  case ASVEC:
    qop_ff_opt.value = 0;
    QOP_asqtad_force_set_opts(&qop_ff_opt, 1);
    eo_fermion_force_asqtad_block( eps, residues, xxx, nterms, veclength );
    break;
  default:  /* FNMAT */
    qop_ff_opt.value = 1;
    QOP_asqtad_force_set_opts(&qop_ff_opt, 1);
    eo_fermion_force_asqtad_multi( eps, residues, xxx, nterms );
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
  default:
    return "FNMAT";
  }
  return NULL;
}

