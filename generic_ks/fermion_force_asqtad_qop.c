/******* fermion_force_asqtad_milc_qop.c ****************/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP fermion force routine */
/* 11/2005 C. DeTar and D. Renner */

/*
 * $Log: fermion_force_asqtad_qop.c,v $
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

static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_ks/fermion_force_asqtad_qop.c,v 1.14 2006/10/12 03:45:16 detar Exp $";

void load_links_and_mom_site(QOP_GaugeField **links, QOP_Force **mom,
			     su3_matrix ***rawlinks, su3_matrix ***rawmom)
{

  /* Copy gauge links from site structure to raw and then to QOP format */
  
  *rawlinks = create_raw_G_from_site_links(EVENANDODD);
  if(*rawlinks == NULL)terminate(1);
  
  *links = QOP_create_G_from_raw((Real **)(*rawlinks),QOP_EVENODD);

  /* Copy momentum from site structure to raw and then to QOP format */

  *rawmom = create_raw_F_from_site_mom(EVENANDODD);
  if(*rawmom == NULL)terminate(1);

  *mom = QOP_create_F_from_raw((Real **)(*rawmom),QOP_EVENODD);
  
}

void unload_links_and_mom_site(QOP_GaugeField **links, QOP_Force **mom,
			       su3_matrix ***rawlinks, su3_matrix ***rawmom)
{

  /* Destroy gauge links and QOP links */

  destroy_raw_G (*rawlinks);   *rawlinks = NULL;
  QOP_destroy_G (*links);      *links = NULL;

  /* Copy momentum from QOP format to raw and then to site structure */

  QOP_extract_F_to_raw((Real **)(*rawmom), *mom, QOP_EVENODD);

  unload_raw_F_to_site_mom(*rawmom, EVENANDODD);

  destroy_raw_F (*rawmom);   *rawmom = NULL;
  QOP_destroy_F (*mom);      *mom = NULL;
}

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
  double dtime;

  dtime=-dclock();
#endif

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("eo_fermion_force_oneterm: Error initializing QOP\n");
    terminate(1);
  }

  /* Load gauge links and momentum */
  load_links_and_mom_site( &links, &mom, &rawlinks, &rawmom );

  /* Copy color vector from site structure to raw and then to QOP format */
  rawvecx = create_raw_V_from_site(x_off,EVENANDODD);
  if(rawvecx == NULL)terminate(1);
  vecx = QOP_create_V_from_raw((Real *)rawvecx,QOP_EVENODD);

  /* Load coefficients */
  load_qop_asqtad_coeffs(&coeff, weight);

  /* Compute fermion force */
  QOP_asqtad_force(&info, links, mom, &coeff, eps, vecx);

  /* Unload momentum and destroy storage for momentum and links */
  unload_links_and_mom_site(  &links, &mom, &rawlinks, &rawmom );

  /* Free QOP source vector */
  QOP_destroy_V(vecx);  vecx = NULL;

  /* Free raw source vector */
  destroy_raw_V(rawvecx); rawvecx = NULL;

#ifdef FFTIME
  dtime += dclock();
node0_printf("FFTIME:  time = %e mflops = %e\n",dtime,
	     (Real)nflop*volume/(1e6*dtime*numnodes()) );
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
  double dtime;

  dtime=-dclock();
#endif

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
  load_qop_asqtad_coeffs(&coeff, 1.);

  /* Compute fermion force */
  QOP_asqtad_force_multi(&info, links, mom, &coeff, epsv, vecx, 2);

  /* Unload momentum and destroy storage for momentum and links */
  unload_links_and_mom_site(  &links, &mom, &rawlinks, &rawmom );

  /* Free QOP source vectors */
  QOP_destroy_V(vecx[0]);  vecx[0] = NULL;
  QOP_destroy_V(vecx[1]);  vecx[1] = NULL;

  /* Free raw source vectors */
  destroy_raw_V(rawvecx[0]); rawvecx[0] = NULL;
  destroy_raw_V(rawvecx[1]); rawvecx[1] = NULL;

#ifdef FFTIME
  dtime += dclock();
  node0_printf("FFTIME:  time = %e mflops = %e\n",dtime,
	       (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif
}


/* Set optimization choice */
/* (For the moment we have only one choice, namely Asqtad!) */ 
enum ks_multiff_opt_t { KS_MULTIFF_ASVEC, KS_MULTIFF_FNMATREV, KS_MULTIFF_FNMAT };
static enum ks_multiff_opt_t ks_multiff_opt = KS_MULTIFF_FNMAT;   /* Default */

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

#ifdef FFTIME
  int nflop = 433968;
  double dtime;

  dtime=-dclock();
#endif

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("eo_fermion_force_multi: Error initializing QOP\n");
    terminate(1);
  }

  /* Load gauge links and momentum */
  load_links_and_mom_site( &links, &mom, &rawlinks, &rawmom );

  /* Make space for QOP vector pointers */
  vecx = (QOP_ColorVector **)malloc(sizeof(QOP_ColorVector *)*nterms);
  if(vecx == NULL){
    printf("eo_fermion_force_multi: no room for vecx pointers\n");
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
  load_qop_asqtad_coeffs(&coeff, 1.0);

  /* Compute fermion force */

  QOP_asqtad_force_multi(&info, links, mom, &coeff, epsv, vecx, nterms);

  /* Unload momentum and destroy storage for momentum and links */
  unload_links_and_mom_site(  &links, &mom, &rawlinks, &rawmom );

  /* Free QOP source vectors */
  for(i = 0; i < nterms; i++){
    QOP_destroy_V(vecx[i]);  vecx[i] = NULL;
  }

  free(epsv);

#ifdef FFTIME
  dtime += dclock();
  node0_printf("FFTIME:  time = %e mflops = %e\n",dtime,
	       (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif

}
