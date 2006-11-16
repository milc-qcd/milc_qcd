/******* fermion_force_asqtad_milc_qop.c ****************/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP fermion force routine */
/* 11/2005 C. DeTar and D. Renner */

/*
 * $Log: fermion_force_asqtad_qop.c,v $
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
#include <qop.h>
#include <string.h>

static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_ks/fermion_force_asqtad_qop.c,v 1.17 2006/11/16 04:17:53 detar Exp $";

/* Standard MILC interface for the single-species Asqtad fermion force routine */
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
  load_qop_asqtad_coeffs(&coeff, weight, get_quark_path_coeff());

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
  node0_printf("FFTIME:  time = %e mflops = %e\n",dtime,
	       (Real)nflop*volume/(1e6*dtime*numnodes()) );
  node0_printf("FFREMAP:  time = %e\n",remaptime);
#endif
}

/* Standard MILC interface for the two-species Asqtad fermion force routine */
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
  load_qop_asqtad_coeffs(&coeff, 1., get_quark_path_coeff());

  /* Compute fermion force */
  dtime = -dclock();
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
  node0_printf("FFTIME:  time = %e mflops = %e\n",dtime,
	       (Real)nflop*volume/(1e6*dtime*numnodes()) );
  node0_printf("FFREMAP:  time = %e\n",remaptime);
#endif
}


/* Set optimization choice */
/* (For the moment we have only one choice, namely Asqtad!) */ 
enum ks_multiff_opt_t { KS_MULTIFF_ASVEC, KS_MULTIFF_FNMATREV, 
			KS_MULTIFF_FNMAT };
static enum ks_multiff_opt_t ks_multiff_opt = KS_MULTIFF_ASVEC;  /* Default */

/* returns 1 for error and 0 for success */

int eo_fermion_force_set_opt(char opt_string[]){
  if(strcmp(opt_string,"ASVEC") == 0)
    ks_multiff_opt = KS_MULTIFF_ASVEC;
  else{
    printf("eo_fermion_force_set_opt: Unrecognized type %s\n",opt_string);
    printf("Currently QOP supports only the ASVEC option\n");
    return 1;
  }
  /*printf("eo_fermion_force set ks_multiff_opt to %d\n",ks_multiff_opt);*/
  return 0;
}


/* Standard MILC interface for multiple terms in the action */

void eo_fermion_force_multi( Real eps, Real *residues, 
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
  char myname[] = "eo_fermion_force_multi";

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
  load_qop_asqtad_coeffs(&coeff, 1., get_quark_path_coeff());

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
  node0_printf("FFTIME:  time = %e mflops = %e\n",dtime,
	       (Real)nflop*volume/(1e6*dtime*numnodes()) );
  node0_printf("FFREMAP:  time = %e\n",remaptime);
#endif

}
