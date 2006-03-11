/******* fermion_force_asqtad_milc_qop.c ****************/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP fermion force routine */
/* 11/2005 C. DeTar and D. Renner */

/*
 * $Log: fermion_force_asqtad_qop.c,v $
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

static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_ks/fermion_force_asqtad_qop.c,v 1.9 2006/03/11 04:24:22 detar Exp $";

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

void load_qop_asqtad_coeffs(QOP_asqtad_coeffs_t *c, int nflavors)
{
  Real *act_path_coeff;
  Real ferm_epsilon;

  /* Load path coefficients from table */
  act_path_coeff = get_quark_path_coeff();

  ferm_epsilon = 2.0*(nflavors/4.0);
  
  /* Path coefficients times fermion epsilon */

  c->one_link     = act_path_coeff[0]*ferm_epsilon ; 
  c->naik         = act_path_coeff[1]*ferm_epsilon ;
  c->three_staple = act_path_coeff[2]*ferm_epsilon ;
  c->five_staple  = act_path_coeff[3]*ferm_epsilon ;
  c->seven_staple = act_path_coeff[4]*ferm_epsilon ;
  c->lepage       = act_path_coeff[5]*ferm_epsilon ;
}

/* Generic MILC interface for the single-species Asqtad fermion force routine */
void eo_fermion_force( Real eps, int nflavors, field_offset x_off )
{

  su3_matrix **rawlinks;
  su3_matrix **rawmom;
  su3_vector *rawvecx;

  QOP_GaugeField *links;
  QOP_Force *mom;
  QOP_ColorVector *vecx;
  
  QOP_asqtad_coeffs_t coeff;

#ifdef FFTIME
  int nflop = 253935;
  double dtime;

  dtime=-dclock();
#endif

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("eo_fermion_force: Error initializing QOP\n");
    terminate(1);
  }

  /* Load gauge links and momentum */
  load_links_and_mom_site( &links, &mom, &rawlinks, &rawmom );

  /* Copy color vector from site structure to raw and then to QOP format */
  rawvecx = create_raw_V_from_site(x_off,EVENANDODD);
  if(rawvecx == NULL)terminate(1);
  vecx = QOP_create_V_from_raw((Real *)rawvecx,QOP_EVENODD);

  /* Load coefficients */
  load_qop_asqtad_coeffs(&coeff, nflavors);

  /* Compute fermion force */
  QOP_asqtad_force(links, mom, &coeff, eps, vecx);

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

/* Generic MILC interface for the two-species Asqtad fermion force routine */
void eo_fermion_force_3f( Real eps, int nflav1, field_offset x1_off, 
			  int nflav2, field_offset x2_off ) {

  su3_matrix **rawlinks;
  su3_matrix **rawmom;
  su3_vector *rawvecx[2];

  QOP_GaugeField *links;
  QOP_Force *mom;
  QOP_ColorVector *vecx[2];
  
  QOP_asqtad_coeffs_t coeff;
  int i;
  Real epsv[2];

#ifdef FFTIME
  int nflop = 433968;
  double dtime;

  dtime=-dclock();
#endif

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("eo_fermion_force: Error initializing QOP\n");
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
  epsv[0] = eps*nflav1/4.0;  epsv[1] = eps*nflav2/4.0;
  load_qop_asqtad_coeffs(&coeff, 4.0);

  /* Compute fermion force */
  QOP_asqtad_force_multi(links, mom, &coeff, epsv, vecx, 2);

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

