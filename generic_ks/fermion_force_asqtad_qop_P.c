/******* fermion_force_asqtad_qop_P.c ****************/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP fermion force routine */
/* 11/2005 C. DeTar and D. Renner */
/* 5/09/07 C. DeTar Support selected precision */
/* NOTE: This is an include file for fermion_force_asqtad_qop_F.c and
   fermion_force_asqtad_qop_D.c.  Any changes must take this usage
   into account */

/* External entry points in this file are private to
   fermion_force_asqtad_qop.c

   eo_fermion_force_oneterm_[FD]
   eo_fermion_force_twoterms_[FD]
   fermion_force_multi_[FD]
   fermion_force_block_[FD]

*/

/*
 * $Log: fermion_force_asqtad_qop_P.c,v $
 * Revision 1.5  2012/11/24 00:02:50  detar
 * Add placeholders for HYPISQ action.  Support HISQ action within ks_imp_dyn.
 *
 * Revision 1.4  2012/02/16 16:30:29  detar
 * Initialize QOP_info_t structure.
 *
 * Revision 1.3  2011/11/29 20:45:56  detar
 * Support new fermion links scheme
 *
 * Revision 1.2  2007/12/14 04:36:31  detar
 * Major modification to support HISQ.
 *
 * Revision 1.1  2007/05/21 05:06:20  detar
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

#if (QOP_PrecisionInt==1)

#define CREATE_F_FROM_SITE4        create_F_F_from_site4
#define CREATE_G_FROM_SITE4        create_F_G_from_site4
#define CREATE_V_FROM_SITE         create_F_V_from_site
#define CREATE_V_FROM_FIELD        create_F_V_from_field
#define EO_FERMION_FORCE_ONETERM   eo_fermion_force_oneterm_F
#define EO_FERMION_FORCE_TWOTERMS  eo_fermion_force_twoterms_F
#define FERMION_FORCE_MULTI        fermion_force_multi_F
#define FERMION_FORCE_BLOCK        fermion_force_block_F
#define MY_REAL                    float
#define LOAD_QOP_ASQTAD_COEFFS     load_qop_F_asqtad_coeffs
#define UNLOAD_F_TO_SITE4          unload_F_F_to_site4
#define UNLOAD_G_TO_SITE4          unload_F_G_to_site4
				   
#else				   
				   
#define CREATE_F_FROM_SITE4        create_D_F_from_site4
#define CREATE_G_FROM_SITE4        create_D_G_from_site4
#define CREATE_V_FROM_SITE         create_D_V_from_site
#define CREATE_V_FROM_FIELD        create_D_V_from_field
#define EO_FERMION_FORCE_ONETERM   eo_fermion_force_oneterm_D
#define EO_FERMION_FORCE_TWOTERMS  eo_fermion_force_twoterms_D
#define FERMION_FORCE_MULTI        fermion_force_multi_D
#define FERMION_FORCE_BLOCK        fermion_force_block_D
#define MY_REAL                    double
#define LOAD_QOP_ASQTAD_COEFFS     load_qop_D_asqtad_coeffs
#define UNLOAD_F_TO_SITE4          unload_D_F_to_site4
#define UNLOAD_G_TO_SITE4          unload_D_G_to_site4
#endif

#include "generic_ks_includes.h"
#include "../include/generic_qop.h"
#include "../include/generic_ks_qop.h"

/* Set default if undeclared */
#ifndef KS_MULTIFF
#define KS_MULTIFF FNMAT
#endif

#ifdef FFTIME
static const char *qop_prec[2] = {"F", "D"};
#endif

//static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_ks/fermion_force_asqtad_qop_P.c,v 1.5 2012/11/24 00:02:50 detar Exp $";

/**********************************************************************/
/* Standard MILC interface for the single-species Asqtad fermion force
   routine */
/**********************************************************************/
void EO_FERMION_FORCE_ONETERM( Real eps, Real weight, su3_vector *x_off,
			       fermion_links_t *fl)
{

  /* For example weight = nflavors/4 */

  QOP_GaugeField *links;
  QOP_Force *mom;
  QOP_ColorVector *vecx;
  
  QOP_asqtad_coeffs_t *coeff = get_action_coeffs(fl);
  QOP_info_t info = {0., 0., 0, 0, 0};

  double remaptime = -dclock();

  QOP_opt_t qop_ff_opt[2] = {
    {.tag = "fnmat_src_min",.value=4},
    {.tag = "veclength",.value=4}
  };

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    node0_printf("eo_fermion_force_oneterm: Error initializing QOP\n");
    terminate(1);
  }

  /* Convert gauge links and momentum to QOP format */

  links = CREATE_G_FROM_SITE4(F_OFFSET(link), EVENANDODD);
  mom   = CREATE_F_FROM_SITE4(F_OFFSET(mom), EVENANDODD);

  /* Convert color vector to QOP format */
  //  vecx = CREATE_V_FROM_SITE(x_off,EVENANDODD);
  vecx = CREATE_V_FROM_FIELD(x_off,EVENANDODD);

//  /* Load coefficients */
//  LOAD_QOP_ASQTAD_COEFFS(&coeff, weight, ap->act_path_coeff);

  /* Compute fermion force */
  remaptime += dclock();
  QOP_asqtad_force_set_opts(qop_ff_opt, 2);
  /* The coefficients are already loaded with weight = 0.5 */
#ifdef QOP_HAS_VERSION // detect newer versions with changed convention
  QOP_asqtad_force(&info, links, mom, coeff, weight*eps, vecx);
#else
  QOP_asqtad_force(&info, links, mom, coeff, 2.*weight*eps, vecx);
#endif
  remaptime -= dclock();

  /* Unload momentum */
  UNLOAD_F_TO_SITE4( F_OFFSET(mom), mom, EVENANDODD );

  /* Clean up */
  QOP_destroy_G(links); links = NULL;
  QOP_destroy_F(mom);   mom = NULL;
  QOP_destroy_V(vecx);  vecx = NULL;

  remaptime += dclock();
#ifdef FFTIME
  node0_printf("FFTIME:  time = %e (qop) terms = 1 mflops = %e\n",
	       info.final_sec, (Real)info.final_flop/(1e6*info.final_sec) );
#ifdef REMAP
  node0_printf("FFREMAP:  time = %e\n",remaptime);
#endif
#endif
}

/**********************************************************************/
/* Standard MILC interface for the two-species Asqtad fermion force
   routine */
/**********************************************************************/
void EO_FERMION_FORCE_TWOTERMS( Real eps, Real weight1, Real weight2, 
				su3_vector *x1_off, su3_vector *x2_off,
				fermion_links_t *fl )
{

  /* For example weight1 = nflavor1/4; weight2 = nflavor2/4 */

  QOP_asqtad_coeffs_t *coeff = get_action_coeffs(fl);
  QOP_GaugeField *links;
  QOP_Force *mom;
  QOP_ColorVector *vecx[2];
  
  MY_REAL epsv[2];
  QOP_info_t info = {0., 0., 0, 0, 0};
  QOP_opt_t qop_ff_opt[2] = {
    {.tag = "fnmat_src_min",.value=4},
    {.tag = "veclength",.value=4}
  };

  double remaptime = -dclock();

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("eo_fermion_force_twoterms: Error initializing QOP\n");
    terminate(1);
  }

  /* Convert gauge links and momentum to QOP format */

  links = CREATE_G_FROM_SITE4(F_OFFSET(link), EVENANDODD);
  mom   = CREATE_F_FROM_SITE4(F_OFFSET(mom), EVENANDODD);

  /* Loop over fermion species */
  /* Convert color vectors in site structure to QOP format */
//  vecx[0] = CREATE_V_FROM_SITE(x1_off,EVENANDODD);
//  vecx[1] = CREATE_V_FROM_SITE(x2_off,EVENANDODD);
  vecx[0] = CREATE_V_FROM_FIELD(x1_off,EVENANDODD);
  vecx[1] = CREATE_V_FROM_FIELD(x2_off,EVENANDODD);

  /* Load coefficients */
#ifdef QOP_HAS_VERSION // detect newer versions with changed convention
  epsv[0] =    eps*weight1;  epsv[1] =    eps*weight2;
#else
  epsv[0] = 2.*eps*weight1;  epsv[1] = 2.*eps*weight2;
#endif
  //  LOAD_QOP_ASQTAD_COEFFS(&coeff, 1., ap->act_path_coeff);

  /* Compute fermion force */
  remaptime += dclock();
  QOP_asqtad_force_set_opts(qop_ff_opt, 2);
  QOP_asqtad_force_multi(&info, links, mom, coeff, epsv, vecx, 2);
  remaptime -= dclock();

  /* Unload momentum */
  UNLOAD_F_TO_SITE4( F_OFFSET(mom), mom, EVENANDODD );

  /* Clean up */
  QOP_destroy_G(links); links = NULL;
  QOP_destroy_F(mom);   mom = NULL;
  QOP_destroy_V(vecx[0]);  vecx[0] = NULL;
  QOP_destroy_V(vecx[1]);  vecx[1] = NULL;

  remaptime += dclock();
#ifdef FFTIME
  node0_printf("FFTIME:  time = %e (qop %s) terms = %d mflops = %e\n",
	       info.final_sec, qop_prec[QOP_PrecisionInt-1],
	       2, (Real)info.final_flop/(1e6*info.final_sec) );
#ifdef REMAP
  node0_printf("FFREMAP:  time = %e\n",remaptime);
#endif
#endif
}


/**********************************************************************/
/*   Version for asqtad.  Parallel transport nterms source vectors    */
/**********************************************************************/

void FERMION_FORCE_MULTI( Real eps, Real *residues, 
			  su3_vector **xxx, int nterms,
			  fermion_links_t *fl ) 
{

  QOP_asqtad_coeffs_t *coeff = get_action_coeffs(fl);
  QOP_GaugeField *links;
  QOP_Force *mom;
  QOP_ColorVector **vecx;
  
  int i;
  MY_REAL *epsv;
  QOP_info_t info = {0., 0., 0, 0, 0};
  char myname[] = "fermion_force_multi";

  double remaptime = -dclock();

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("%s(%d): Error initializing QOP\n",myname,this_node);
    terminate(1);
  }

  /* Convert gauge links and momentum to QOP format */

  links = CREATE_G_FROM_SITE4(F_OFFSET(link), EVENANDODD);
  mom   = CREATE_F_FROM_SITE4(F_OFFSET(mom), EVENANDODD);

  /* Make space for QOP vector pointers */
  vecx = (QOP_ColorVector **)malloc(sizeof(QOP_ColorVector *)*nterms);
  if(vecx == NULL){
    printf("%s(%d): no room for vecx pointers\n",myname,this_node);
    terminate(1);
  }

  /* Convert color vectors to QOP format */
  for(i = 0; i < nterms; i++){
    vecx[i] = CREATE_V_FROM_FIELD(xxx[i], EVENANDODD);
  }

  /* Make space for weights */
  epsv = (MY_REAL *)malloc(sizeof(MY_REAL)*nterms);
  /* Load coefficients */
#ifdef QOP_HAS_VERSION // detect newer versions with changed convention
  for(i = 0; i < nterms; i++) epsv[i] = eps*residues[i];
#else
  for(i = 0; i < nterms; i++) epsv[i] = 2.*eps*residues[i];
#endif
  //  LOAD_QOP_ASQTAD_COEFFS(&coeff, 1., ap->act_path_coeff);

  /* Compute fermion force */
  remaptime += dclock();
  QOP_asqtad_force_multi(&info, links, mom, coeff, epsv, vecx, nterms);
  remaptime -= dclock();

  /* Unload momentum */
  UNLOAD_F_TO_SITE4( F_OFFSET(mom), mom, EVENANDODD );

  /* Clean up */
  QOP_destroy_G(links); links = NULL;
  QOP_destroy_F(mom);   mom = NULL;
  for(i = 0; i < nterms; i++){
    QOP_destroy_V(vecx[i]);  vecx[i] = NULL;
  }
  free(vecx);  vecx = NULL;
  free(epsv);  epsv = NULL;

  remaptime += dclock();
#ifdef FFTIME
  node0_printf("FFTIME:  time = %e (qop %s) terms = %d mflops = %e\n",
	       info.final_sec, qop_prec[QOP_PrecisionInt-1],
	       nterms, info.final_flop/(1e6*info.final_sec) );
#ifdef REMAP
  node0_printf("FFREMAP:  time = %e\n",remaptime);
#endif
#endif

}


/**********************************************************************/
/*   Parallel transport vectors in blocks of veclength. Asqtad only   */
/**********************************************************************/

void FERMION_FORCE_BLOCK( Real eps, Real *residues, 
			  su3_vector **xxx, int nterms, int veclength,
			  fermion_links_t *fl) {

  int j;

  /* First do blocks of size veclength */
  for( j = 0;  j <= nterms-veclength; j += veclength )
    FERMION_FORCE_MULTI( eps, &(residues[j]), xxx+j, veclength, fl );
  
#ifndef ONEMASS
  /* Continue with pairs if needed */
  for( ; j <= nterms-2 ; j+=2 ){
    EO_FERMION_FORCE_TWOTERMS( eps, residues[j], residues[j+1],
			       xxx[j], xxx[j + 1], fl );
  }

  /* Finish with a single if needed */
  for( ; j <= nterms-1; j++ ){
    EO_FERMION_FORCE_ONETERM( eps, residues[j], xxx[j], fl );
  }
#else
  /* Thrown in just to make it compile.  Probably never used. */
  /* Finish with a single if needed */
  for( ; j <= nterms-1; j++ ){
    EO_FERMION_FORCE_ONETERM( eps, residues[j], xxx[j], fl );
  }
#endif

}

